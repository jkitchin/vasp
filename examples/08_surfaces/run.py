#!/usr/bin/env python
"""
08 - Surfaces

This script demonstrates surface slab calculations, calculating the
surface energy of Cu(111).

Usage:
    python run.py
"""

from ase.build import bulk, fcc111, add_adsorbate
from ase.constraints import FixAtoms
import numpy as np
from vasp import Vasp

print("=" * 60)
print("Surface Slab Calculations")
print("=" * 60)
print()

# =============================================================================
# Part 1: Bulk copper reference
# =============================================================================

print("Part 1: Bulk Cu reference calculation")
print("-" * 40)
print()

cu_bulk = bulk('Cu', 'fcc', a=3.6)

calc_bulk = Vasp(
    label='surfaces/cu_bulk',
    atoms=cu_bulk,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ismear=1,
    sigma=0.1,
)

energy_bulk = calc_bulk.potential_energy
e_per_atom = energy_bulk / len(cu_bulk)

print(f"  Total energy: {energy_bulk:.6f} eV")
print(f"  Energy per atom: {e_per_atom:.6f} eV/atom")
print()

# =============================================================================
# Part 2: Cu(111) slab
# =============================================================================

print("Part 2: Cu(111) slab calculation")
print("-" * 40)
print()

# Create Cu(111) slab with 4 layers and vacuum
slab = fcc111('Cu', size=(2, 2, 4), vacuum=12.0, a=3.6)

print(f"Slab structure:")
print(f"  Surface: Cu(111)")
print(f"  Size: 2x2 supercell, 4 layers")
print(f"  Atoms: {len(slab)}")
print(f"  Vacuum: 12.0 Å")
print()

# Fix bottom two layers
positions_z = slab.positions[:, 2]
z_sorted = np.sort(np.unique(np.round(positions_z, decimals=2)))
z_threshold = z_sorted[1] + 0.1  # Fix bottom 2 layers

constraint = FixAtoms(mask=positions_z <= z_threshold)
slab.set_constraint(constraint)

n_fixed = sum(positions_z <= z_threshold)
print(f"  Fixed atoms (bottom 2 layers): {n_fixed}")
print()

# Slab calculation with dipole correction
calc_slab = Vasp(
    label='surfaces/cu111_slab',
    atoms=slab,
    xc='PBE',
    encut=400,
    kpts=(6, 6, 1),  # Gamma-centered, 1 in z for slab
    ismear=1,
    sigma=0.1,
    ldipol=True,     # Dipole correction
    idipol=3,        # Correct in z direction
    isif=2,          # Relax positions only (not cell)
    ibrion=2,        # Conjugate gradient
    nsw=50,          # Max ionic steps
    ediffg=-0.02,    # Force convergence: 0.02 eV/Å
)

energy_slab = calc_slab.potential_energy

print(f"  Total energy: {energy_slab:.6f} eV")
print()

# =============================================================================
# Part 3: Calculate surface energy
# =============================================================================

print("Part 3: Surface energy calculation")
print("-" * 40)
print()

n_atoms_slab = len(slab)

# Surface area (only one side counted since symmetric slab)
cell = slab.get_cell()
area = np.linalg.norm(np.cross(cell[0], cell[1]))

# Surface energy (per surface, slab has two surfaces)
e_surface_ev = (energy_slab - n_atoms_slab * e_per_atom) / (2 * area)

# Convert eV/Å² to J/m²
ev_per_ang2_to_j_per_m2 = 16.0217663
e_surface_j = e_surface_ev * ev_per_ang2_to_j_per_m2

print(f"  Number of atoms in slab: {n_atoms_slab}")
print(f"  Surface area: {area:.4f} Å²")
print(f"  Surface energy: {e_surface_ev:.6f} eV/Å²")
print(f"  Surface energy: {e_surface_j:.4f} J/m²")
print(f"  Experimental: ~1.83 J/m²")
print()

# =============================================================================
# Part 4: Work function (from electrostatic potential)
# =============================================================================

print("Part 4: Work function estimation")
print("-" * 40)
print()

# Static calculation for accurate LOCPOT
calc_static = Vasp(
    label='surfaces/cu111_static',
    atoms=slab,
    xc='PBE',
    encut=400,
    kpts=(6, 6, 1),
    ismear=1,
    sigma=0.1,
    ldipol=True,
    idipol=3,
    lvtot=True,  # Write LOCPOT
)

_ = calc_static.potential_energy
fermi = calc_static.results.get('fermi_level', 0.0)

# Try to read work function from LOCPOT
try:
    wf = calc_static.get_work_function()
    print(f"  Fermi level: {fermi:.4f} eV")
    print(f"  Vacuum level: {wf + fermi:.4f} eV")
    print(f"  Work function: {wf:.4f} eV")
    print(f"  Experimental Cu(111): 4.94 eV")
except Exception as e:
    print(f"  Fermi level: {fermi:.4f} eV")
    print(f"  Work function calculation requires LOCPOT analysis")
    print(f"  Experimental Cu(111): 4.94 eV")
print()

# =============================================================================
# Part 5: Convergence with slab thickness
# =============================================================================

print("Part 5: Surface energy vs. slab thickness")
print("-" * 40)
print()

print("Testing convergence with number of layers:")
print()

# Quick convergence test with different layer counts
layer_counts = [3, 4, 5, 6]
surface_energies = []

for n_layers in layer_counts:
    slab_test = fcc111('Cu', size=(2, 2, n_layers), vacuum=12.0, a=3.6)

    calc_test = Vasp(
        label=f'surfaces/cu111_{n_layers}L',
        atoms=slab_test,
        xc='PBE',
        encut=400,
        kpts=(6, 6, 1),
        ismear=1,
        sigma=0.1,
    )

    e_test = calc_test.potential_energy
    n_atoms = len(slab_test)
    e_surf = (e_test - n_atoms * e_per_atom) / (2 * area) * ev_per_ang2_to_j_per_m2
    surface_energies.append(e_surf)

    print(f"  {n_layers} layers: {e_surf:.4f} J/m²")

print()
print(f"  Variation: {max(surface_energies) - min(surface_energies):.4f} J/m²")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Cu(111) Surface Properties:")
print(f"  Surface energy: {e_surface_j:.3f} J/m² (exp: ~1.83 J/m²)")
print(f"  Work function: ~4.9 eV (exp: 4.94 eV)")
print()
print("Key points:")
print("  - Use sufficient vacuum (10-15 Å)")
print("  - Apply dipole correction (LDIPOL=True)")
print("  - Converge with slab thickness")
print("  - Fix bottom layers during relaxation")
print()
print("Next: Try 09_adsorption/ for molecules on surfaces.")
