#!/usr/bin/env python
"""
14 - Van der Waals Corrections

This script demonstrates dispersion-corrected DFT calculations using
graphite as an example, where van der Waals interactions are essential.

Standard DFT fails to describe London dispersion forces, leading to
incorrect binding in layered materials, molecular crystals, and
adsorption systems.

Usage:
    python run.py
"""

from ase import Atoms
from ase.build import bulk
import numpy as np
from vasp import Vasp
from vasp.parameters import get_vdw_params

print("=" * 60)
print("Van der Waals Corrections")
print("=" * 60)
print()

# =============================================================================
# Part 1: Graphite structure
# =============================================================================

print("Part 1: Graphite structure")
print("-" * 40)
print()

# Create graphite (AB stacking)
# Hexagonal cell with 4 atoms
a = 2.46   # In-plane lattice constant
c = 6.71   # Interlayer spacing (experimental)

graphite = Atoms(
    symbols=['C', 'C', 'C', 'C'],
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [1/3, 2/3, 0.0],
        [0.0, 0.0, 0.5],
        [2/3, 1/3, 0.5],
    ],
    cell=[
        [a, 0, 0],
        [-a/2, a*np.sqrt(3)/2, 0],
        [0, 0, c],
    ],
    pbc=True,
)

print("Structure: Graphite (AB stacking)")
print(f"  In-plane lattice: a = {a} Å")
print(f"  Interlayer spacing: c/2 = {c/2:.2f} Å")
print(f"  Experimental interlayer: 3.35 Å")
print()

# =============================================================================
# Part 2: Standard PBE (no vdW)
# =============================================================================

print("Part 2: Standard PBE (no dispersion)")
print("-" * 40)
print()

calc_pbe = Vasp(
    label='vdw/graphite_pbe',
    atoms=graphite.copy(),
    xc='PBE',
    encut=500,
    kpts=(12, 12, 6),
    ismear=0,
    sigma=0.05,
)

e_pbe = calc_pbe.potential_energy
print(f"PBE energy: {e_pbe:.6f} eV")
print()

# Scan interlayer distance with PBE
print("Interlayer binding scan (PBE):")
c_values = [5.5, 6.0, 6.5, 7.0, 8.0, 10.0]
pbe_energies = []

for c_test in c_values:
    graphite_test = graphite.copy()
    cell = graphite_test.get_cell()
    cell[2, 2] = c_test
    graphite_test.set_cell(cell, scale_atoms=True)

    calc = Vasp(
        label=f'vdw/graphite_pbe_c{c_test:.1f}',
        atoms=graphite_test,
        xc='PBE',
        encut=500,
        kpts=(12, 12, 6),
        ismear=0,
        sigma=0.05,
    )
    e = calc.potential_energy
    pbe_energies.append(e)
    print(f"  c = {c_test:.1f} Å: E = {e:.4f} eV")

print()
print("PBE shows no binding minimum (energy decreases with c)")
print("This is WRONG - graphite layers should bind!")
print()

# =============================================================================
# Part 3: D3-BJ correction
# =============================================================================

print("Part 3: PBE-D3(BJ) (Grimme's D3 with BJ damping)")
print("-" * 40)
print()

d3bj_params = get_vdw_params('d3bj')
print("D3-BJ parameters:")
print(f"  IVDW = 12 (D3 with BJ damping)")
print(f"  Empirical C6 coefficients with geometry-dependent terms")
print()

d3bj_energies = []
for c_test in c_values:
    graphite_test = graphite.copy()
    cell = graphite_test.get_cell()
    cell[2, 2] = c_test
    graphite_test.set_cell(cell, scale_atoms=True)

    calc = Vasp(
        label=f'vdw/graphite_d3bj_c{c_test:.1f}',
        atoms=graphite_test,
        xc='PBE',
        encut=500,
        kpts=(12, 12, 6),
        ismear=0,
        sigma=0.05,
        **d3bj_params,
    )
    e = calc.potential_energy
    d3bj_energies.append(e)
    print(f"  c = {c_test:.1f} Å: E = {e:.4f} eV")

# Find minimum
min_idx = np.argmin(d3bj_energies)
print()
print(f"D3-BJ minimum at c ≈ {c_values[min_idx]:.1f} Å")
print(f"Experimental: c = 6.71 Å (interlayer = 3.35 Å)")
print()

# =============================================================================
# Part 4: Compare dispersion methods
# =============================================================================

print("Part 4: Comparison of vdW methods")
print("-" * 40)
print()

# Test at experimental c
graphite_exp = graphite.copy()

methods = [
    ('None', {}),
    ('d2', get_vdw_params('d2')),
    ('d3', get_vdw_params('d3')),
    ('d3bj', get_vdw_params('d3bj')),
    ('ts', get_vdw_params('ts')),
]

print("Comparison at experimental geometry (c = 6.71 Å):")
print()
print(f"  {'Method':<10} {'Energy (eV)':<15} {'IVDW':<10}")
print("  " + "-" * 35)

for name, params in methods:
    calc = Vasp(
        label=f'vdw/graphite_{name}',
        atoms=graphite_exp.copy(),
        xc='PBE',
        encut=500,
        kpts=(12, 12, 6),
        ismear=0,
        sigma=0.05,
        **params,
    )
    e = calc.potential_energy
    ivdw = params.get('ivdw', 0)
    print(f"  {name:<10} {e:<15.4f} {ivdw:<10}")

print()

# =============================================================================
# Part 5: vdW-DF functionals
# =============================================================================

print("Part 5: vdW-DF nonlocal functionals")
print("-" * 40)
print()

print("VASP supports nonlocal vdW density functionals:")
print()
print("  vdW-DF (LUSE_VDW=True, GGA=RE):")
print("    - Nonlocal correlation from electron density")
print("    - More physical than empirical D2/D3")
print("    - More expensive computationally")
print()
print("  vdW-DF2, optB86b-vdW, optPBE-vdW, etc.")
print()

# These require LUSE_VDW=True and specific GGA settings
# Not demonstrated here due to computational cost

# =============================================================================
# Part 6: Binding energy calculation
# =============================================================================

print("Part 6: Graphite binding energy")
print("-" * 40)
print()

# Calculate isolated graphene layer
graphene = Atoms(
    symbols=['C', 'C'],
    scaled_positions=[
        [0.0, 0.0, 0.5],
        [1/3, 2/3, 0.5],
    ],
    cell=[
        [a, 0, 0],
        [-a/2, a*np.sqrt(3)/2, 0],
        [0, 0, 15.0],  # Large vacuum
    ],
    pbc=True,
)

calc_graphene = Vasp(
    label='vdw/graphene_d3bj',
    atoms=graphene,
    xc='PBE',
    encut=500,
    kpts=(12, 12, 1),
    ismear=0,
    sigma=0.05,
    **d3bj_params,
)

e_graphene = calc_graphene.potential_energy

# Binding energy per atom
e_binding = (e_pbe - 2 * e_graphene) / 4  # Per atom
e_binding_d3 = (d3bj_energies[c_values.index(6.5)] - 2 * e_graphene) / 4

print("Interlayer binding energy (per C atom):")
print()
print(f"  PBE (no vdW): ~0 meV (unbound)")
print(f"  PBE-D3(BJ): ~{e_binding_d3*1000:.1f} meV")
print(f"  Experimental: ~52 meV")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()

print("Available dispersion corrections:")
print()
print("  Empirical (pairwise):")
print("    - D2: Original Grimme, fixed C6")
print("    - D3: Geometry-dependent C6")
print("    - D3-BJ: D3 with Becke-Johnson damping (recommended)")
print("    - TS: Tkatchenko-Scheffler, uses electron density")
print("    - TS-SCS: TS with self-consistent screening")
print()
print("  Nonlocal (density-based):")
print("    - vdW-DF family: Physical but expensive")
print("    - SCAN+rVV10: Meta-GGA with nonlocal correlation")
print()

print("Recommendations:")
print("  - D3-BJ: Good accuracy, low cost (default choice)")
print("  - TS-SCS: Better for metals and surfaces")
print("  - vdW-DF: When electron density matters")
print()
print("Key points:")
print("  - Standard DFT misses dispersion completely")
print("  - Essential for: layered materials, molecules, adsorption")
print("  - Different methods have different accuracy/cost tradeoffs")
print()
print("Next: Try 15_workflows/ for automated calculation pipelines.")
