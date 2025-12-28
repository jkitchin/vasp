#!/usr/bin/env python
"""
19 - Pseudopotential Selection

This script demonstrates how to choose and specify VASP pseudopotentials
(POTCARs) for different types of calculations.

Usage:
    python run.py
"""

from ase import Atoms
from ase.build import bulk

from vasp import Vasp

# =============================================================================
# Overview
# =============================================================================

print("=" * 60)
print("VASP Pseudopotential Selection Guide")
print("=" * 60)
print()
print("This example demonstrates how to choose the right POTCAR")
print("for your VASP calculations using the 'setups' parameter.")
print()

# =============================================================================
# Example 1: Standard vs _pv for Transition Metals
# =============================================================================

print("-" * 60)
print("Example 1: Transition Metal (Fe)")
print("-" * 60)
print()

# Create BCC Fe
fe = bulk('Fe', 'bcc', a=2.87)
fe.set_initial_magnetic_moments([2.5])

print("Standard Fe POTCAR (8 valence electrons):")
print("  - Treats 3d and 4s as valence")
print("  - Freezes 3s and 3p core electrons")
print()

# Standard POTCAR - no setups needed
calc_standard = Vasp(
    label='fe_standard',
    atoms=fe.copy(),
    xc='PBE',
    encut=400,
    kpts=(12, 12, 12),
    ispin=2,
    # No setups = use default Fe POTCAR
)

print("Fe_pv POTCAR (14 valence electrons):")
print("  - Includes 3p semi-core states")
print("  - More accurate for oxides and magnetic systems")
print("  - Requires higher ENCUT (~450 eV)")
print()

# Fe_pv POTCAR
calc_pv = Vasp(
    label='fe_pv',
    atoms=fe.copy(),
    xc='PBE',
    encut=450,  # Higher cutoff for _pv
    kpts=(12, 12, 12),
    ispin=2,
    setups={'Fe': 'pv'},  # Use Fe_pv POTCAR
)

print("Python code:")
print("  # Standard:")
print("  calc = Vasp('fe', atoms=fe, encut=400)")
print()
print("  # With _pv:")
print("  calc = Vasp('fe', atoms=fe, encut=450, setups={'Fe': 'pv'})")
print()

# =============================================================================
# Example 2: Perovskite Oxides (ABO3)
# =============================================================================

print("-" * 60)
print("Example 2: Perovskite Oxide (SrTiO3)")
print("-" * 60)
print()

# Create cubic SrTiO3
a = 3.905
srtio3 = Atoms(
    'SrTiO3',
    positions=[
        (0.0, 0.0, 0.0),      # Sr
        (a/2, a/2, a/2),      # Ti
        (a/2, a/2, 0.0),      # O
        (a/2, 0.0, a/2),      # O
        (0.0, a/2, a/2),      # O
    ],
    cell=[a, a, a],
    pbc=True,
)

print("Recommended setups for perovskites:")
print("  Sr: _sv (includes 4s4p, 10 valence e-)")
print("  Ti: _pv (includes 3p, 10 valence e-)")
print("  O:  standard (6 valence e-)")
print()

calc_perov = Vasp(
    label='srtio3',
    atoms=srtio3,
    xc='PBE',
    encut=520,  # Higher cutoff for _sv/_pv
    kpts=(6, 6, 6),
    setups={
        'Sr': 'sv',  # Include semi-core
        'Ti': 'pv',  # Include 3p states
        'O': '',     # Standard is fine
    },
)

print("Python code:")
print("  setups = {'Sr': 'sv', 'Ti': 'pv', 'O': ''}")
print("  calc = Vasp('srtio3', atoms=srtio3, encut=520, setups=setups)")
print()

# =============================================================================
# Example 3: Battery Materials
# =============================================================================

print("-" * 60)
print("Example 3: Battery Material (LiFePO4)")
print("-" * 60)
print()

print("Recommended setups for Li-ion cathodes:")
print("  Li: _sv (CRITICAL - includes 1s, 3 valence e-)")
print("  Fe: _pv (includes 3p, 14 valence e-)")
print("  P:  standard (5 valence e-)")
print("  O:  standard (6 valence e-)")
print()
print("Why Li_sv is important:")
print("  - Standard Li has only 1 valence electron")
print("  - Li_sv treats 1s as valence (3 electrons)")
print("  - Critical for accurate Li intercalation voltages")
print()

print("Python code:")
print("  setups = {'Li': 'sv', 'Fe': 'pv', 'P': '', 'O': ''}")
print("  calc = Vasp('lifepo4', atoms=atoms, encut=520, setups=setups)")
print()

# =============================================================================
# Example 4: ENCUT Recommendations
# =============================================================================

print("-" * 60)
print("Example 4: ENCUT Recommendations")
print("-" * 60)
print()

print("Each POTCAR has an ENMAX (recommended minimum cutoff).")
print("Rule: ENCUT >= 1.3 × max(ENMAX of all elements)")
print()

encut_table = """
Element   Standard  _pv    _sv    Recommended ENCUT
---------------------------------------------------------
Li        140       -      499    650 eV (with _sv)
Na        102       265    646    520 eV (with _pv)
Fe        268       293    391    520 eV (with _pv)
Ti        178       222    275    520 eV (with _pv)
O         400       -      -      520 eV
Sr        -         -      229    520 eV
"""
print(encut_table)

print("Typical settings:")
print("  Standard elements only:    ENCUT = 520 eV")
print("  With _pv variants:         ENCUT = 520 eV")
print("  With _sv (Li, Na):         ENCUT = 650+ eV")
print()

# =============================================================================
# Summary: Quick Reference
# =============================================================================

print("=" * 60)
print("Quick Reference")
print("=" * 60)
print()

print("POTCAR Variants:")
print("  (none)  Standard - default, computationally efficient")
print("  _pv     p-valence - includes semi-core p states")
print("  _sv     Semi-core valence - includes more core states")
print("  _d      d-electrons - includes d for main group")
print("  _h      Hard - higher cutoff, more accurate")
print("  _s      Soft - lower cutoff, faster")
print("  _GW     GW-ready - for GW calculations")
print()

print("When to use _pv:")
print("  - 3d transition metals in oxides")
print("  - Magnetic systems")
print("  - High-pressure calculations")
print()

print("When to use _sv:")
print("  - Alkali metals (Li, Na, K)")
print("  - Alkaline earth metals (Ca, Sr, Ba)")
print("  - Battery and energy storage materials")
print()

print("Specifying setups in code:")
print("  from vasp import Vasp")
print()
print("  calc = Vasp(")
print("      'my_calc',")
print("      atoms=atoms,")
print("      setups={")
print("          'Fe': 'pv',   # Use Fe_pv")
print("          'Li': 'sv',   # Use Li_sv")
print("          'O': '',      # Use standard O")
print("      },")
print("      encut=520,")
print("  )")
print()

print("=" * 60)
print("Done! See README.md for complete documentation.")
print("=" * 60)
print()
print("Next: Review the examples to understand the trade-offs")
print("      between accuracy and computational cost.")
