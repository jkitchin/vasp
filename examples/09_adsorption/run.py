#!/usr/bin/env python
"""
09 - Adsorption

This script demonstrates calculating adsorption energies for molecules
on metal surfaces, using CO on Pt(111) as an example.

Usage:
    python run.py
"""

from ase.build import bulk, fcc111, add_adsorbate, molecule
from ase.constraints import FixAtoms
import numpy as np
from vasp import Vasp

print("=" * 60)
print("Adsorption Energy Calculations")
print("=" * 60)
print()

# =============================================================================
# Part 1: Reference calculations
# =============================================================================

print("Part 1: Reference calculations")
print("-" * 40)
print()

# Gas-phase CO molecule
print("1a. Gas-phase CO molecule")
co = molecule('CO')
co.center(vacuum=8.0)

calc_co = Vasp(
    label='adsorption/co_gas',
    atoms=co,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),  # Molecule in box - gamma only
    ismear=0,
    sigma=0.05,
    ispin=2,         # CO is spin-polarized in gas phase
)

e_co = calc_co.potential_energy
print(f"  CO energy: {e_co:.6f} eV")
print()

# Pt bulk reference
print("1b. Pt bulk reference")
pt_bulk = bulk('Pt', 'fcc', a=3.92)

calc_pt_bulk = Vasp(
    label='adsorption/pt_bulk',
    atoms=pt_bulk,
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ismear=1,
    sigma=0.1,
)

e_pt_bulk = calc_pt_bulk.potential_energy
e_pt_per_atom = e_pt_bulk / len(pt_bulk)
print(f"  Pt bulk energy per atom: {e_pt_per_atom:.6f} eV")
print()

# =============================================================================
# Part 2: Clean Pt(111) surface
# =============================================================================

print("Part 2: Clean Pt(111) surface")
print("-" * 40)
print()

# Create Pt(111) slab
# 2x2 supercell, 4 layers, 12 Å vacuum
slab = fcc111('Pt', size=(2, 2, 4), vacuum=12.0, a=3.92)

# Fix bottom 2 layers
positions_z = slab.positions[:, 2]
z_sorted = np.sort(np.unique(np.round(positions_z, decimals=2)))
z_threshold = z_sorted[1] + 0.1
constraint = FixAtoms(mask=positions_z <= z_threshold)
slab.set_constraint(constraint)

print(f"  Slab: Pt(111) 2x2, 4 layers")
print(f"  Atoms: {len(slab)}")
print()

calc_slab = Vasp(
    label='adsorption/pt111_clean',
    atoms=slab,
    xc='PBE',
    encut=400,
    kpts=(4, 4, 1),
    ismear=1,
    sigma=0.1,
    ldipol=True,
    idipol=3,
    isif=2,
    ibrion=2,
    nsw=30,
    ediffg=-0.03,
)

e_slab = calc_slab.potential_energy
print(f"  Clean slab energy: {e_slab:.6f} eV")
print()

# =============================================================================
# Part 3: CO adsorption at different sites
# =============================================================================

print("Part 3: CO adsorption at different sites")
print("-" * 40)
print()

adsorption_sites = ['ontop', 'bridge', 'fcc', 'hcp']
adsorption_results = {}

for site in adsorption_sites:
    print(f"  Calculating {site} site...")

    # Create fresh slab with CO
    slab_co = fcc111('Pt', size=(2, 2, 4), vacuum=12.0, a=3.92)

    # Fix bottom layers
    positions_z = slab_co.positions[:, 2]
    z_sorted = np.sort(np.unique(np.round(positions_z, decimals=2)))
    z_threshold = z_sorted[1] + 0.1
    constraint = FixAtoms(mask=positions_z <= z_threshold)
    slab_co.set_constraint(constraint)

    # Add CO at specified site
    # height is approximate initial guess
    add_adsorbate(slab_co, 'CO', height=1.9, position=site)

    calc_ads = Vasp(
        label=f'adsorption/pt111_co_{site}',
        atoms=slab_co,
        xc='PBE',
        encut=400,
        kpts=(4, 4, 1),
        ismear=1,
        sigma=0.1,
        ldipol=True,
        idipol=3,
        isif=2,
        ibrion=2,
        nsw=50,
        ediffg=-0.03,
    )

    e_ads = calc_ads.potential_energy

    # Calculate adsorption energy
    e_adsorption = e_ads - e_slab - e_co

    adsorption_results[site] = {
        'energy': e_ads,
        'e_ads': e_adsorption,
    }

    print(f"    Total energy: {e_ads:.6f} eV")
    print(f"    Adsorption energy: {e_adsorption:.3f} eV")
    print()

# =============================================================================
# Part 4: Analysis
# =============================================================================

print("Part 4: Summary of adsorption energies")
print("-" * 40)
print()

print(f"{'Site':<10} {'E_ads (eV)':<12} {'Stability':<15}")
print("-" * 37)

# Sort by adsorption energy (more negative = more stable)
sorted_sites = sorted(adsorption_results.items(), key=lambda x: x[1]['e_ads'])

for i, (site, data) in enumerate(sorted_sites):
    stability = "Most stable" if i == 0 else f"ΔE = {data['e_ads'] - sorted_sites[0][1]['e_ads']:.3f} eV"
    print(f"{site:<10} {data['e_ads']:<12.3f} {stability:<15}")

print()

# =============================================================================
# Part 5: Coverage effects (conceptual)
# =============================================================================

print("Part 5: Coverage dependence")
print("-" * 40)
print()

print("Different coverages can be studied by changing supercell size:")
print()
print(f"  {'Coverage':<15} {'Supercell':<12} {'CO per cell':<12}")
print("  " + "-" * 39)
print(f"  {'1 ML':<15} {'1x1':<12} {'1':<12}")
print(f"  {'0.5 ML':<15} {'2x1':<12} {'1':<12}")
print(f"  {'0.25 ML':<15} {'2x2':<12} {'1':<12}")
print(f"  {'0.11 ML':<15} {'3x3':<12} {'1':<12}")
print()

print("Coverage effects include:")
print("  - Lateral adsorbate-adsorbate interactions")
print("  - Substrate-mediated interactions")
print("  - Work function changes")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary: CO/Pt(111) Adsorption")
print("=" * 60)
print()

most_stable = sorted_sites[0]
print(f"Most stable site: {most_stable[0]}")
print(f"Adsorption energy: {most_stable[1]['e_ads']:.3f} eV")
print()
print("Experimental data:")
print("  - CO prefers atop site on Pt(111)")
print("  - Experimental E_ads: ~-1.4 to -1.5 eV")
print("  - Note: Standard DFT may overestimate binding")
print()
print("Key points:")
print("  - E_ads = E(slab+adsorbate) - E(slab) - E(adsorbate)")
print("  - Negative values indicate favorable adsorption")
print("  - Test different sites to find global minimum")
print("  - Consider coverage effects for real systems")
print()
print("Next: Try 10_reactions/ for reaction energetics.")
