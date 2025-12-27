#!/usr/bin/env python
"""
01 - Getting Started: Your First VASP Calculation

This script calculates the total energy of bulk silicon.
It's the simplest possible VASP calculation to verify your setup works.

Usage:
    python run.py
"""

from ase.build import bulk
from vasp import Vasp

# =============================================================================
# Step 1: Create the atomic structure
# =============================================================================

# Silicon in the diamond structure
# The experimental lattice constant is 5.43 Å
atoms = bulk('Si', 'diamond', a=5.43)

print("Structure created:")
print(f"  Formula: {atoms.get_chemical_formula()}")
print(f"  Number of atoms: {len(atoms)}")
print(f"  Lattice constant: {atoms.cell[0, 1]:.3f} Å")
print()

# =============================================================================
# Step 2: Set up the VASP calculator
# =============================================================================

calc = Vasp(
    # Directory for calculation files
    label='si_energy',

    # Attach the atoms
    atoms=atoms,

    # Exchange-correlation functional
    xc='PBE',

    # Plane-wave energy cutoff (eV)
    # Start low for testing, increase for production
    encut=300,

    # k-point mesh
    # Denser grids give more accurate results but take longer
    kpts=(4, 4, 4),

    # Smearing method and width
    # ismear=1 (Methfessel-Paxton) is good for metals/semiconductors
    # sigma should be small enough that entropy term is < 1 meV/atom
    ismear=1,
    sigma=0.1,

    # Don't write large files for this simple test
    lwave=False,   # Don't write WAVECAR
    lcharg=False,  # Don't write CHGCAR
)

# =============================================================================
# Step 3: Run the calculation
# =============================================================================

print("Running VASP calculation...")
print("=" * 50)

# This triggers the actual VASP run
energy = calc.potential_energy

print("=" * 50)
print()

# =============================================================================
# Step 4: Print results
# =============================================================================

print("Results:")
print(f"  Total energy: {energy:.6f} eV")
print(f"  Energy per atom: {energy / len(atoms):.6f} eV/atom")
print()

# Access additional results
fermi = calc.results.get('fermi_level')
if fermi:
    print(f"  Fermi level: {fermi:.4f} eV")

# Check convergence
if calc.results.get('converged', True):
    print("  Status: Converged successfully")
else:
    print("  Status: WARNING - Not converged!")

print()
print(f"Output files written to: {calc.directory}/")
print()
print("Congratulations! Your first VASP calculation is complete.")
print("Next: Try 02_convergence/ to learn about parameter convergence.")
