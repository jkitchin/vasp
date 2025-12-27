#!/usr/bin/env python
"""
10 - Reactions

This script demonstrates calculating reaction energetics:
- Dissociation energies
- Reaction pathways
- Transition state concepts

Using H2 dissociation as a simple example.

Usage:
    python run.py
"""

from ase.build import molecule
from ase import Atoms
import numpy as np
from vasp import Vasp

print("=" * 60)
print("Reaction Energetics")
print("=" * 60)
print()

# =============================================================================
# Part 1: Reference energies
# =============================================================================

print("Part 1: Reference molecule energies")
print("-" * 40)
print()

# H2 molecule
print("1a. H2 molecule")
h2 = molecule('H2')
h2.center(vacuum=8.0)

calc_h2 = Vasp(
    label='reactions/h2',
    atoms=h2,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ismear=0,
    sigma=0.05,
    ispin=2,
)

e_h2 = calc_h2.potential_energy
print(f"  H2 energy: {e_h2:.6f} eV")
print()

# Single H atom (for dissociation reference)
print("1b. H atom")
h_atom = Atoms('H', positions=[[0, 0, 0]])
h_atom.center(vacuum=8.0)

calc_h = Vasp(
    label='reactions/h_atom',
    atoms=h_atom,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[1.0],  # H atom is spin-1/2
)

e_h = calc_h.potential_energy
print(f"  H atom energy: {e_h:.6f} eV")
print()

# Dissociation energy
e_diss = 2 * e_h - e_h2
print(f"  H2 dissociation energy: {e_diss:.3f} eV")
print(f"  Experimental: 4.52 eV")
print()

# =============================================================================
# Part 2: O2 and H2O for combustion example
# =============================================================================

print("Part 2: Combustion reaction (2H2 + O2 → 2H2O)")
print("-" * 40)
print()

# O2 molecule (triplet ground state)
print("2a. O2 molecule")
o2 = molecule('O2')
o2.center(vacuum=8.0)

calc_o2 = Vasp(
    label='reactions/o2',
    atoms=o2,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[1.0, 1.0],  # O2 is triplet
)

e_o2 = calc_o2.potential_energy
print(f"  O2 energy: {e_o2:.6f} eV")
print()

# H2O molecule
print("2b. H2O molecule")
h2o = molecule('H2O')
h2o.center(vacuum=8.0)

calc_h2o = Vasp(
    label='reactions/h2o',
    atoms=h2o,
    xc='PBE',
    encut=400,
    kpts=(1, 1, 1),
    ismear=0,
    sigma=0.05,
    ispin=2,
)

e_h2o = calc_h2o.potential_energy
print(f"  H2O energy: {e_h2o:.6f} eV")
print()

# Combustion reaction: 2H2 + O2 -> 2H2O
e_reaction = 2 * e_h2o - 2 * e_h2 - e_o2
print("Reaction: 2H2 + O2 → 2H2O")
print(f"  Reaction energy: {e_reaction:.3f} eV")
print(f"  Per H2O formed: {e_reaction/2:.3f} eV")
print(f"  Experimental: ~-2.5 eV per H2O")
print()

# =============================================================================
# Part 3: Potential energy surface scan
# =============================================================================

print("Part 3: H2 dissociation curve")
print("-" * 40)
print()

print("Scanning H-H bond distance...")
print()

bond_lengths = [0.6, 0.74, 0.9, 1.1, 1.4, 1.8, 2.5, 3.5]
energies = []

for d in bond_lengths:
    h2_stretched = Atoms('H2',
        positions=[[0, 0, 0], [d, 0, 0]])
    h2_stretched.center(vacuum=8.0)

    calc_scan = Vasp(
        label=f'reactions/h2_d{d:.2f}',
        atoms=h2_stretched,
        xc='PBE',
        encut=400,
        kpts=(1, 1, 1),
        ismear=0,
        sigma=0.05,
        ispin=2,
    )

    e = calc_scan.potential_energy
    energies.append(e)
    print(f"  d = {d:.2f} Å: E = {e:.4f} eV")

print()

# Find equilibrium
min_idx = np.argmin(energies)
print(f"  Minimum at d = {bond_lengths[min_idx]:.2f} Å")
print(f"  Experimental: 0.74 Å")
print()

# =============================================================================
# Part 4: Transition state concepts
# =============================================================================

print("Part 4: Transition state concepts")
print("-" * 40)
print()

print("Finding transition states requires specialized methods:")
print()
print("  1. NEB (Nudged Elastic Band):")
print("     - Connect reactant and product states")
print("     - Find minimum energy path")
print("     - Identify saddle point (TS)")
print()
print("  2. Dimer Method:")
print("     - Start near suspected TS")
print("     - Climb to nearest saddle point")
print()
print("  3. CI-NEB (Climbing Image NEB):")
print("     - Enhanced NEB with better TS convergence")
print()
print("VASP parameters for NEB:")
print("  IMAGES = N       (number of images)")
print("  SPRING = -5      (spring constant)")
print("  IBRION = 1,3     (optimizer)")
print("  LCLIMB = .TRUE.  (climbing image)")
print()

# =============================================================================
# Part 5: Activation energy example (conceptual)
# =============================================================================

print("Part 5: Reaction energy diagram")
print("-" * 40)
print()

# Conceptual example: A + B -> AB
# (Using H + H -> H2 energies)

e_reactants = 2 * e_h  # Two separate H atoms
e_product = e_h2       # H2 molecule
e_barrier = e_reactants + 0.1  # Hypothetical small barrier

print("Example: H + H → H2")
print()
print(f"  Reactants (2H): {e_reactants:.3f} eV")
print(f"  TS (barrier):   {e_barrier:.3f} eV")
print(f"  Product (H2):   {e_product:.3f} eV")
print()
print(f"  Activation energy (Ea): {e_barrier - e_reactants:.3f} eV")
print(f"  Reaction energy (ΔE):   {e_product - e_reactants:.3f} eV")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Key reaction energetics:")
print(f"  H2 dissociation: {e_diss:.2f} eV (exp: 4.52 eV)")
print(f"  H2 combustion: {e_reaction/2:.2f} eV/H2O (exp: ~-2.5 eV)")
print()
print("Important concepts:")
print("  - Reaction energy = E(products) - E(reactants)")
print("  - Activation energy = E(TS) - E(reactants)")
print("  - Zero-point corrections can be significant")
print("  - Entropic contributions for gas-phase")
print()
print("Methods for finding transition states:")
print("  - NEB: Robust but computationally expensive")
print("  - Dimer: Fast but needs good initial guess")
print("  - CI-NEB: Best for accurate barrier heights")
print()
print("Next: Try 11_phonons/ for vibrational properties.")
