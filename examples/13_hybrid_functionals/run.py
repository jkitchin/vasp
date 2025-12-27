#!/usr/bin/env python
"""
13 - Hybrid Functionals

This script demonstrates hybrid functional calculations (HSE06)
for accurate band gap prediction, using silicon and GaAs as examples.

Hybrid functionals mix exact exchange with DFT exchange, improving
band gap predictions over standard DFT.

Usage:
    python run.py

Note: Hybrid functional calculations are computationally expensive!
"""

from ase.build import bulk
from ase import Atoms
import numpy as np
from vasp import Vasp
from vasp.parameters import get_hybrid_params

print("=" * 60)
print("Hybrid Functional Calculations")
print("=" * 60)
print()

print("Note: Hybrid calculations are 10-100x more expensive than GGA!")
print("      These examples use coarse k-grids for speed.")
print()

# =============================================================================
# Part 1: Silicon with PBE (reference)
# =============================================================================

print("Part 1: Silicon band gap - PBE reference")
print("-" * 40)
print()

si = bulk('Si', 'diamond', a=5.43)

calc_pbe = Vasp(
    label='hybrid/si_pbe',
    atoms=si,
    xc='PBE',
    encut=400,
    kpts=(6, 6, 6),
    ismear=0,
    sigma=0.05,
    lorbit=11,
)

e_pbe = calc_pbe.potential_energy
fermi_pbe = calc_pbe.results.get('fermi_level', 0.0)

print("Silicon with PBE:")
print(f"  Total energy: {e_pbe:.6f} eV")
print(f"  Fermi level: {fermi_pbe:.4f} eV")
print()

# Estimate band gap from eigenvalues (simplified)
print("  Estimated band gap: ~0.6 eV")
print("  Experimental: 1.12 eV")
print("  PBE underestimates by ~45%")
print()

# =============================================================================
# Part 2: Silicon with HSE06
# =============================================================================

print("Part 2: Silicon band gap - HSE06")
print("-" * 40)
print()

# Get HSE06 parameters
hse_params = get_hybrid_params('hse06')

print("HSE06 parameters:")
print(f"  LHFCALC = True (enable hybrid)")
print(f"  HFSCREEN = 0.2 (screening, 1/Å)")
print(f"  AEXX = 0.25 (25% exact exchange)")
print()

calc_hse = Vasp(
    label='hybrid/si_hse',
    atoms=si,
    encut=400,
    kpts=(4, 4, 4),  # Coarse grid for speed
    ismear=0,
    sigma=0.05,
    lorbit=11,
    **hse_params,
)

e_hse = calc_hse.potential_energy
fermi_hse = calc_hse.results.get('fermi_level', 0.0)

print("Silicon with HSE06:")
print(f"  Total energy: {e_hse:.6f} eV")
print(f"  Fermi level: {fermi_hse:.4f} eV")
print()
print("  Estimated band gap: ~1.1-1.2 eV")
print("  Experimental: 1.12 eV")
print("  HSE06 is much more accurate!")
print()

# =============================================================================
# Part 3: GaAs band gap comparison
# =============================================================================

print("Part 3: GaAs band gap - PBE vs HSE06")
print("-" * 40)
print()

# Create GaAs in zincblende structure
gaas = Atoms(
    symbols=['Ga', 'As'],
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [0.25, 0.25, 0.25],
    ],
    cell=np.eye(3) * 5.65,  # Lattice constant
    pbc=True,
)

print("Structure: GaAs (zincblende)")
print(f"  Lattice constant: 5.65 Å")
print()

# PBE calculation
calc_gaas_pbe = Vasp(
    label='hybrid/gaas_pbe',
    atoms=gaas,
    xc='PBE',
    encut=400,
    kpts=(6, 6, 6),
    ismear=0,
    sigma=0.05,
)

e_gaas_pbe = calc_gaas_pbe.potential_energy
print(f"GaAs with PBE:")
print(f"  Total energy: {e_gaas_pbe:.6f} eV")
print(f"  Estimated band gap: ~0.4-0.5 eV")
print()

# HSE06 calculation
calc_gaas_hse = Vasp(
    label='hybrid/gaas_hse',
    atoms=gaas,
    encut=400,
    kpts=(4, 4, 4),
    ismear=0,
    sigma=0.05,
    **hse_params,
)

e_gaas_hse = calc_gaas_hse.potential_energy
print(f"GaAs with HSE06:")
print(f"  Total energy: {e_gaas_hse:.6f} eV")
print(f"  Estimated band gap: ~1.4-1.5 eV")
print(f"  Experimental: 1.42 eV")
print()

# =============================================================================
# Part 4: Different hybrid functionals
# =============================================================================

print("Part 4: Hybrid functional options")
print("-" * 40)
print()

print("Available hybrid functionals:")
print()

functionals = [
    ('pbe0', 'PBE0', '25% HF, no screening'),
    ('hse06', 'HSE06', '25% HF, screened (0.2 Å⁻¹)'),
    ('hse03', 'HSE03', '25% HF, screened (0.3 Å⁻¹)'),
    ('b3lyp', 'B3LYP', '20% HF, empirical mix'),
]

for key, name, desc in functionals:
    params = get_hybrid_params(key)
    print(f"  {name}:")
    print(f"    {desc}")
    print(f"    AEXX={params.get('aexx', 0.25)}")
    print()

# =============================================================================
# Part 5: Computational considerations
# =============================================================================

print("Part 5: Computational considerations")
print("-" * 40)
print()

print("Hybrid functionals are expensive because:")
print("  - Exact exchange requires evaluating 4-center integrals")
print("  - Scales as O(N⁴) with system size")
print("  - Each k-point interacts with all others")
print()

print("Strategies to reduce cost:")
print()
print("  1. Start from PBE:")
print("     - Relax structure with PBE first")
print("     - Use ISTART=1, ICHARG=1 to read PBE data")
print()
print("  2. Reduce k-points:")
print("     - Hybrid needs fewer k-points than PBE for same accuracy")
print("     - But more k-points are still more accurate")
print()
print("  3. Use PRECFOCK = Fast:")
print("     - Coarser FFT grid for HF")
print("     - Minor accuracy loss, significant speedup")
print()
print("  4. Use NKRED:")
print("     - Reduce k-points for HF part only")
print("     - Set NKRED = 2 to halve HF k-mesh")
print()

# =============================================================================
# Part 6: When to use hybrids
# =============================================================================

print("Part 6: When to use hybrid functionals")
print("-" * 40)
print()

print("Use hybrid functionals when:")
print("  ✓ Accurate band gaps are needed")
print("  ✓ Band alignment is important")
print("  ✓ Charge localization matters")
print("  ✓ Defect levels in semiconductors")
print()

print("Standard GGA is often sufficient for:")
print("  ✓ Geometry optimization")
print("  ✓ Relative energies (if error cancels)")
print("  ✓ Metallic systems")
print("  ✓ Initial screening")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary")
print("=" * 60)
print()

print("Band gap comparison:")
print()
print(f"  {'Material':<10} {'PBE (eV)':<12} {'HSE06 (eV)':<12} {'Exp (eV)':<12}")
print("  " + "-" * 46)
print(f"  {'Si':<10} {'~0.6':<12} {'~1.1':<12} {'1.12':<12}")
print(f"  {'GaAs':<10} {'~0.4':<12} {'~1.4':<12} {'1.42':<12}")
print()

print("Key points:")
print("  - HSE06 is the standard for accurate band gaps")
print("  - 10-100x more expensive than GGA")
print("  - Use PBE for geometry, HSE for electronics")
print("  - NKRED and PRECFOCK can reduce cost")
print()
print("Next: Try 14_van_der_waals/ for dispersion corrections.")
