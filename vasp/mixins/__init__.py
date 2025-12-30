"""Mixin classes providing VASP calculator functionality.

These mixins are combined with the base calculator to provide
the full VASP interface. Each mixin groups related methods:

- IOMixin: File reading and writing (INCAR, POSCAR, KPOINTS, etc.)
- ElectronicMixin: Electronic structure properties
- AnalysisMixin: Post-processing analysis (Bader, band structure)
- DynamicsMixin: NEB and vibrational calculations
"""

from .analysis import AnalysisMixin
from .dynamics import DynamicsMixin
from .electronic import ElectronicMixin
from .io import IOMixin

__all__ = [
    'IOMixin',
    'ElectronicMixin',
    'AnalysisMixin',
    'DynamicsMixin',
]
