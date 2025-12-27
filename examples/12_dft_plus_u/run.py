#!/usr/bin/env python
"""
12 - DFT+U

This script demonstrates DFT+U calculations for strongly correlated
systems, using iron oxide (FeO, wustite) as an example.

DFT+U adds a Hubbard U correction to better describe localized d or f
electrons that standard DFT handles poorly.

Usage:
    python run.py
"""

from ase import Atoms
from ase.build import bulk
import numpy as np
from vasp import Vasp
from vasp.parameters import get_ldau_params, HubbardU

print("=" * 60)
print("DFT+U Calculations")
print("=" * 60)
print()

# =============================================================================
# Part 1: Standard DFT for FeO
# =============================================================================

print("Part 1: Standard DFT calculation")
print("-" * 40)
print()

# Create FeO in rock salt structure
a_FeO = 4.33  # Experimental lattice constant
feo = Atoms(
    symbols=['Fe', 'O'],
    scaled_positions=[
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
    ],
    cell=[a_FeO, a_FeO, a_FeO],
    pbc=True,
)

print("Structure: Rock salt FeO (Wustite)")
print(f"  Lattice constant: {a_FeO} Å")
print(f"  Fe: 3d⁶ configuration")
print()

# Standard DFT calculation (spin-polarized)
calc_dft = Vasp(
    label='dft_u/feo_standard',
    atoms=feo,
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[4.0, 0.0],  # Fe high-spin, O non-magnetic
    lorbit=11,
)

e_dft = calc_dft.potential_energy
magmom_dft = calc_dft.results.get('magnetic_moment', 0.0)
magmoms_dft = calc_dft.results.get('magnetic_moments', [])

print("Standard DFT results:")
print(f"  Total energy: {e_dft:.6f} eV")
print(f"  Total magnetic moment: {magmom_dft:.4f} μB")
if len(magmoms_dft) > 0:
    print(f"  Fe magnetic moment: {magmoms_dft[0]:.4f} μB")
print()

# =============================================================================
# Part 2: DFT+U for FeO
# =============================================================================

print("Part 2: DFT+U calculation")
print("-" * 40)
print()

# Get DFT+U parameters for Fe-O system
# Typical U values for Fe in oxides: 4-5 eV
ldau_params = get_ldau_params(
    symbols=['Fe', 'O'],
    u_values={'Fe': HubbardU(u=4.0, j=0.0)},
)

print("DFT+U parameters:")
print(f"  U(Fe d-electrons) = 4.0 eV")
print(f"  J(Fe) = 0.0 eV (Dudarev formulation: U_eff = U - J)")
print()

calc_dftu = Vasp(
    label='dft_u/feo_u4',
    atoms=feo.copy(),
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ismear=0,
    sigma=0.05,
    ispin=2,
    magmom=[4.0, 0.0],
    lorbit=11,
    **ldau_params,
)

e_dftu = calc_dftu.potential_energy
magmom_dftu = calc_dftu.results.get('magnetic_moment', 0.0)
magmoms_dftu = calc_dftu.results.get('magnetic_moments', [])

print("DFT+U (U=4.0 eV) results:")
print(f"  Total energy: {e_dftu:.6f} eV")
print(f"  Total magnetic moment: {magmom_dftu:.4f} μB")
if len(magmoms_dftu) > 0:
    print(f"  Fe magnetic moment: {magmoms_dftu[0]:.4f} μB")
    print(f"  Experimental: ~3.6 μB")
print()

# =============================================================================
# Part 3: U parameter scan
# =============================================================================

print("Part 3: Effect of U parameter")
print("-" * 40)
print()

print("Scanning U values...")
print()

u_values = [0.0, 2.0, 4.0, 5.0, 6.0]
results = []

for u in u_values:
    ldau = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=u)})

    calc = Vasp(
        label=f'dft_u/feo_u{u:.1f}',
        atoms=feo.copy(),
        xc='PBE',
        encut=500,
        kpts=(8, 8, 8),
        ismear=0,
        sigma=0.05,
        ispin=2,
        magmom=[4.0, 0.0],
        lorbit=11,
        **ldau,
    )

    e = calc.potential_energy
    mag = calc.results.get('magnetic_moments', [0.0])[0]
    results.append({'u': u, 'energy': e, 'magmom': mag})

    print(f"  U = {u:.1f} eV: Fe moment = {mag:.3f} μB")

print()
print("Effect of U on Fe magnetic moment:")
print(f"  {'U (eV)':<10} {'μ_Fe (μB)':<12}")
print("  " + "-" * 22)
for r in results:
    print(f"  {r['u']:<10.1f} {r['magmom']:<12.3f}")
print()

# =============================================================================
# Part 4: Compare band gap
# =============================================================================

print("Part 4: Electronic structure comparison")
print("-" * 40)
print()

print("Band gap comparison:")
print()

# Note: Full band gap analysis requires DOS calculation
# Here we show the concept

print("Standard DFT:")
print("  - FeO predicted to be metallic (incorrect)")
print("  - Underestimated d-electron localization")
print()

print("DFT+U (U=4 eV):")
print("  - Correct insulating behavior")
print("  - Band gap opens due to d-electron localization")
print("  - Experimental gap: ~2.4 eV")
print()

# =============================================================================
# Part 5: Different U schemes
# =============================================================================

print("Part 5: DFT+U formulations")
print("-" * 40)
print()

print("VASP DFT+U options (LDAUTYPE):")
print()
print("  1. Liechtenstein (LDAUTYPE=1):")
print("     - Uses both U and J separately")
print("     - E_U = (U-J)/2 × Σ [n_m - Σ n_mm' n_m'm]")
print()
print("  2. Dudarev (LDAUTYPE=2) [default]:")
print("     - Uses effective U_eff = U - J")
print("     - E_U = (U-J)/2 × Σ [n_m(1 - n_m)]")
print("     - Simpler, most commonly used")
print()
print("  4. Liechtenstein without exchange (LDAUTYPE=4):")
print("     - Uses U without J")
print()

# =============================================================================
# Part 6: How to determine U
# =============================================================================

print("Part 6: Determining U values")
print("-" * 40)
print()

print("Methods to determine Hubbard U:")
print()
print("  1. Linear response (constrained DFT):")
print("     - Self-consistent calculation of U")
print("     - Most rigorous approach")
print()
print("  2. Empirical fitting:")
print("     - Match to experimental properties")
print("     - Band gap, magnetic moment, lattice constant")
print()
print("  3. Literature values:")
print("     - Use established U values for similar systems")
print()
print("Common U values (eV):")
print("  Fe(d): 3-5   Co(d): 3-4   Ni(d): 5-7")
print("  Mn(d): 3-4   Cr(d): 3-4   V(d): 3-4")
print("  Ce(f): 5-6   U(f): 4-5")
print()

# =============================================================================
# Summary
# =============================================================================

print("=" * 60)
print("Summary: DFT+U for FeO")
print("=" * 60)
print()

print("Comparison:")
print(f"  {'Method':<15} {'Fe moment (μB)':<18} {'Behavior':<15}")
print("  " + "-" * 48)
print(f"  {'Standard DFT':<15} {magmoms_dft[0] if magmoms_dft else 0.0:<18.3f} {'Metallic':<15}")
print(f"  {'DFT+U (U=4)':<15} {magmoms_dftu[0] if magmoms_dftu else 0.0:<18.3f} {'Insulating':<15}")
print(f"  {'Experimental':<15} {'~3.6':<18} {'Insulating':<15}")
print()

print("Key points:")
print("  - DFT+U corrects for self-interaction error")
print("  - Essential for transition metal oxides")
print("  - U value affects magnetic moment and band gap")
print("  - Dudarev formulation (LDAUTYPE=2) is most common")
print()
print("Next: Try 13_hybrid_functionals/ for HSE06 calculations.")
