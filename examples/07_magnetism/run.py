#!/usr/bin/env python
"""
07 - Magnetism

This script demonstrates spin-polarized calculations for magnetic materials,
using iron (Fe) and nickel oxide (NiO) as examples.

Usage:
    python run.py
"""

from ase.build import bulk
from ase import Atoms
import numpy as np
from vasp import Vasp

print("=" * 60)
print("Magnetic Calculations")
print("=" * 60)
print()

# =============================================================================
# Part 1: Ferromagnetic Iron (BCC)
# =============================================================================

print("Part 1: Ferromagnetic Iron")
print("-" * 40)
print()

# Create BCC iron
fe = bulk('Fe', 'bcc', a=2.87)

print("Structure: BCC Fe")
print(f"  Atoms: {len(fe)}")
print(f"  Lattice constant: 2.87 Å")
print()

# Ferromagnetic calculation
calc_fe = Vasp(
    label='mag/fe_fm',
    atoms=fe,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ismear=1,           # Methfessel-Paxton for metals
    sigma=0.1,
    ispin=2,            # Spin-polarized
    magmom=[3.0],       # Initial magnetic moment per atom
    lorbit=11,          # For magnetic moment output
)

energy_fe = calc_fe.potential_energy
magmom_fe = calc_fe.results.get('magnetic_moment', 0.0)
magmoms_fe = calc_fe.results.get('magnetic_moments', [])

print("Ferromagnetic Fe results:")
print(f"  Total energy: {energy_fe:.6f} eV")
print(f"  Total magnetic moment: {magmom_fe:.4f} μB")
if len(magmoms_fe) > 0:
    print(f"  Magnetic moment per atom: {magmoms_fe[0]:.4f} μB")
print(f"  Experimental: 2.22 μB/atom")
print()

# =============================================================================
# Part 2: Non-magnetic Fe for comparison
# =============================================================================

print("Part 2: Non-magnetic Fe (for comparison)")
print("-" * 40)
print()

calc_fe_nm = Vasp(
    label='mag/fe_nm',
    atoms=fe,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ismear=1,
    sigma=0.1,
    ispin=1,  # Non-spin-polarized
)

energy_fe_nm = calc_fe_nm.potential_energy

print("Non-magnetic Fe results:")
print(f"  Total energy: {energy_fe_nm:.6f} eV")
print()

delta_e = energy_fe_nm - energy_fe
print(f"Magnetic stabilization energy: {delta_e:.3f} eV")
print("  (Ferromagnetic state is lower in energy)")
print()

# =============================================================================
# Part 3: Antiferromagnetic NiO (rock salt)
# =============================================================================

print("Part 3: Antiferromagnetic NiO")
print("-" * 40)
print()

# Create NiO in rock salt structure
# Need a 2-atom unit cell for antiferromagnetic ordering
a_NiO = 4.17  # Lattice constant
nio = Atoms(
    symbols=['Ni', 'O'],
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ],
    cell=[a_NiO, a_NiO, a_NiO],
    pbc=True,
)

print("Structure: Rock salt NiO")
print(f"  Lattice constant: {a_NiO} Å")
print()

# For AFM ordering in NiO, we need a larger supercell
# The [111] AFM ordering requires doubling along [111]
# Here we use a simple ferromagnetic calculation for demonstration
# (full AFM requires larger supercell)

print("Ferromagnetic NiO (demonstration):")

calc_nio = Vasp(
    label='mag/nio_fm',
    atoms=nio,
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[2.0, 0.0],  # Initial: Ni magnetic, O non-magnetic
    lorbit=11,
)

energy_nio = calc_nio.potential_energy
magmom_nio = calc_nio.results.get('magnetic_moment', 0.0)
magmoms_nio = calc_nio.results.get('magnetic_moments', [])

print(f"  Total energy: {energy_nio:.6f} eV")
print(f"  Total magnetic moment: {magmom_nio:.4f} μB")
if len(magmoms_nio) > 0:
    print(f"  Ni magnetic moment: {magmoms_nio[0]:.4f} μB")
    print(f"  O magnetic moment: {magmoms_nio[1]:.4f} μB")
print()

# =============================================================================
# Part 4: NiO with DFT+U
# =============================================================================

print("Part 4: NiO with DFT+U correction")
print("-" * 40)
print()

from vasp.parameters import get_ldau_params, HubbardU

# Get DFT+U parameters for Ni-O system
ldau_params = get_ldau_params(
    symbols=['Ni', 'O'],
    u_values={'Ni': HubbardU(u=6.45, j=0.0)},  # Standard U for Ni in oxides
)

print(f"DFT+U parameters: U(Ni) = 6.45 eV")
print()

calc_nio_u = Vasp(
    label='mag/nio_u',
    atoms=nio,
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[2.0, 0.0],
    lorbit=11,
    **ldau_params,
)

energy_nio_u = calc_nio_u.potential_energy
magmom_nio_u = calc_nio_u.results.get('magnetic_moment', 0.0)
magmoms_nio_u = calc_nio_u.results.get('magnetic_moments', [])

print(f"  Total energy: {energy_nio_u:.6f} eV")
print(f"  Total magnetic moment: {magmom_nio_u:.4f} μB")
if len(magmoms_nio_u) > 0:
    print(f"  Ni magnetic moment: {magmoms_nio_u[0]:.4f} μB")
    print(f"  Experimental Ni moment: ~1.9 μB")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Fe (BCC):")
print(f"  Magnetic moment: {magmom_fe:.2f} μB (exp: 2.22 μB)")
print(f"  Magnetic stabilization: {delta_e:.3f} eV")
print()
print("NiO (rock salt):")
print(f"  DFT magnetic moment: {magmom_nio:.2f} μB")
print(f"  DFT+U magnetic moment: {magmom_nio_u:.2f} μB")
print()
print("Key points:")
print("  - ISPIN=2 enables spin polarization")
print("  - MAGMOM sets initial magnetic moments")
print("  - DFT+U improves description of correlated d/f electrons")
print("  - AFM ordering requires appropriate supercells")
print()
print("Next: Try 08_surfaces/ for slab calculations.")
