#!/usr/bin/env python
"""
03 - Structure Relaxation

This script demonstrates geometry optimization in VASP,
including both ionic relaxation and full cell relaxation.

Usage:
    python run.py
"""

from ase.build import bulk
from ase import Atoms
import numpy as np
from vasp import Vasp

# =============================================================================
# Create a distorted silicon structure
# =============================================================================

print("=" * 60)
print("Structure Relaxation Example")
print("=" * 60)
print()

# Start with experimental lattice constant (we'll let VASP find the PBE value)
atoms = bulk('Si', 'diamond', a=5.43)

# Slightly distort the structure to give the optimizer something to do
# Move one atom off its ideal position
atoms.positions[1] += [0.1, 0.05, -0.08]

print("Initial structure:")
print(f"  Lattice constant: {atoms.cell[0, 1]:.4f} Å")
print(f"  Atom 1 position: {atoms.positions[0]}")
print(f"  Atom 2 position: {atoms.positions[1]}")
print()

# =============================================================================
# Part 1: Ionic Relaxation (ISIF=2)
# =============================================================================

print("=" * 60)
print("Part 1: Ionic Relaxation (ISIF=2)")
print("=" * 60)
print()
print("Relaxing atomic positions at fixed cell...")
print()

calc_ionic = Vasp(
    label='relax_ionic',
    atoms=atoms.copy(),
    xc='PBE',
    encut=400,
    kpts=(4, 4, 4),
    ismear=1,
    sigma=0.1,

    # Relaxation settings
    ibrion=2,      # Conjugate gradient
    isif=2,        # Relax ions only, not cell
    nsw=50,        # Maximum steps
    ediffg=-0.02,  # Force convergence: |F| < 0.02 eV/Å

    lwave=False,
    lcharg=False,
)

energy_ionic = calc_ionic.potential_energy
forces_ionic = calc_ionic.results.get('forces')
max_force_ionic = np.max(np.abs(forces_ionic)) if forces_ionic is not None else 0

print(f"Results after ionic relaxation:")
print(f"  Energy: {energy_ionic:.6f} eV")
print(f"  Maximum force: {max_force_ionic:.4f} eV/Å")
print(f"  Relaxed positions:")
print(f"    Atom 1: {calc_ionic.atoms.positions[0]}")
print(f"    Atom 2: {calc_ionic.atoms.positions[1]}")
print()

# =============================================================================
# Part 2: Full Cell Relaxation (ISIF=3)
# =============================================================================

print("=" * 60)
print("Part 2: Full Cell Relaxation (ISIF=3)")
print("=" * 60)
print()
print("Relaxing both atoms and cell parameters...")
print()

# Start fresh with a slightly wrong lattice constant
atoms_full = bulk('Si', 'diamond', a=5.50)  # 1% too large

calc_full = Vasp(
    label='relax_full',
    atoms=atoms_full.copy(),
    xc='PBE',
    encut=400,
    kpts=(4, 4, 4),
    ismear=1,
    sigma=0.1,

    # Relaxation settings
    ibrion=2,      # Conjugate gradient
    isif=3,        # Relax ions AND cell (volume + shape)
    nsw=50,        # Maximum steps
    ediffg=-0.02,  # Force convergence

    lwave=False,
    lcharg=False,
)

energy_full = calc_full.potential_energy
forces_full = calc_full.results.get('forces')
max_force_full = np.max(np.abs(forces_full)) if forces_full is not None else 0

# Calculate relaxed lattice constant
# For diamond structure, a = 2 * cell[0,1] / sqrt(2) ... but simpler to just read cell
relaxed_cell = calc_full.atoms.get_cell()
# For FCC-like cell, extract lattice constant
a_relaxed = relaxed_cell[0, 1] * 2  # Approximate

print(f"Results after full relaxation:")
print(f"  Energy: {energy_full:.6f} eV")
print(f"  Maximum force: {max_force_full:.4f} eV/Å")
print(f"  Initial lattice constant: 5.50 Å")
print(f"  Relaxed cell parameter: {relaxed_cell[0, 1]:.4f} Å")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()
print("Relaxation type comparison:")
print(f"  Ionic only (ISIF=2): E = {energy_ionic:.4f} eV")
print(f"  Full cell (ISIF=3):  E = {energy_full:.4f} eV")
print()
print("Key points:")
print("  - Use ISIF=2 for surfaces (fixed vacuum)")
print("  - Use ISIF=3 for bulk structures")
print("  - PBE typically overestimates lattice constants by ~1%")
print("  - Check OUTCAR to verify convergence")
print()
print("Next: Try 04_equation_of_state/ to calculate bulk modulus.")
