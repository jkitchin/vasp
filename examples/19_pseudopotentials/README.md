# Example 19: Pseudopotential Selection

## Overview

Choosing the right pseudopotential (POTCAR) is critical for accurate VASP calculations. This example covers:

- Types of VASP pseudopotentials and their naming conventions
- When to use each variant
- How to specify custom setups in the calculator
- ENMAX considerations for plane-wave cutoff

## VASP Pseudopotential Types

VASP uses PAW (Projector Augmented Wave) pseudopotentials. For each element, multiple versions are available:

### Standard POTCARs

The default POTCAR (e.g., `Fe`, `O`, `Si`) treats only valence electrons explicitly. These are suitable for most calculations and are computationally efficient.

### Variant Suffixes

| Suffix | Meaning | When to Use |
|--------|---------|-------------|
| `_sv` | Semi-core valence | Include next shell as valence (e.g., 3s3p for Fe) |
| `_pv` | p-valence | Include p semi-core states |
| `_d` | d-electrons | Include d states in valence for main group elements |
| `_h` | Hard | Larger cutoff, more accurate, slower |
| `_s` | Soft | Lower cutoff, less accurate, faster |
| `_GW` | GW-ready | For GW calculations, more accurate core |

### Common Examples

| Element | Standard | Common Variants | Notes |
|---------|----------|-----------------|-------|
| Fe | Fe (8 valence) | Fe_pv (14), Fe_sv (16) | Use _pv for accurate magnetism |
| Ti | Ti (4 valence) | Ti_pv (10), Ti_sv (12) | Use _pv for oxides |
| O | O (6 valence) | O_h, O_s | O_h for high accuracy |
| Si | Si (4 valence) | Si_h, Si_d | Standard is usually fine |
| Li | Li (1 valence) | Li_sv (3) | Use _sv for batteries |
| Na | Na (1 valence) | Na_pv (7), Na_sv (9) | Use _pv minimum |
| K | K (9 valence) | K_sv (9) | Standard already includes 3s3p |
| Ca | Ca (2 valence) | Ca_pv (8), Ca_sv (10) | Use _pv or _sv for perovskites |

## Guidelines for Selection

### 1. Transition Metals

**Recommendation**: Use `_pv` variants for 3d transition metals.

The standard potentials freeze the 3s3p electrons, which can cause errors when:
- Computing magnetic properties
- Studying oxides or other ionic compounds
- High-pressure calculations

```python
# Good practice for transition metal oxides
calc = Vasp(
    'fe2o3',
    atoms=atoms,
    setups={'Fe': 'pv'},  # Include 3p states
    encut=520,            # Higher cutoff needed
)
```

### 2. Alkali and Alkaline Earth Metals

**Recommendation**: Use `_pv` or `_sv` variants.

The standard potentials have very few valence electrons, leading to:
- Poor description of bonding
- Errors in lattice constants
- Issues with charge transfer

```python
# For battery materials with Li
calc = Vasp(
    'lifepo4',
    atoms=atoms,
    setups={'Li': 'sv', 'Fe': 'pv'},
    encut=520,
)
```

### 3. Main Group Elements

**Recommendation**: Standard potentials are usually sufficient.

For O, N, C, Si, etc., the standard potentials work well. Use `_h` (hard) variants only when:
- Extremely high accuracy is needed
- Studying core-level properties
- The element is under high pressure

### 4. Rare Earths and Actinides

**Recommendation**: Use `_3` or appropriate f-electron variant.

These elements have complex electronic structures. Consult the VASP wiki for specific recommendations.

## ENMAX Considerations

Each POTCAR has an ENMAX value (recommended minimum ENCUT). Using variant potentials typically requires higher cutoffs:

| Variant | Typical ENMAX Increase |
|---------|------------------------|
| Standard | Reference |
| `_pv` | +50-100 eV |
| `_sv` | +100-200 eV |
| `_h` | +100-150 eV |
| `_GW` | +50-100 eV |

**Best Practice**: Set ENCUT to 1.3× the maximum ENMAX across all elements:

```python
# The calculator automatically reads ENMAX from POTCARs
# You can set ENCUT explicitly for reproducibility
calc = Vasp(
    'calculation',
    atoms=atoms,
    setups={'Ti': 'pv', 'O': 'h'},
    encut=600,  # Higher than default due to variants
)
```

## Specifying Setups in the Calculator

### Method 1: Dictionary

```python
from vasp import Vasp

calc = Vasp(
    'srtio3',
    atoms=atoms,
    setups={
        'Sr': 'sv',   # Uses Sr_sv POTCAR
        'Ti': 'pv',   # Uses Ti_pv POTCAR
        'O': '',      # Uses standard O POTCAR (explicit)
    },
)
```

### Method 2: Recommended Presets

For common material classes, use parameter presets:

```python
from vasp.parameters import get_recommended_setups

# Get recommended setups for a list of elements
setups = get_recommended_setups(['Fe', 'O'], accuracy='high')
# Returns: {'Fe': 'pv', 'O': ''}

calc = Vasp('fe2o3', atoms=atoms, setups=setups)
```

## Material-Specific Recommendations

### Perovskite Oxides (ABO3)

```python
setups = {
    'Sr': 'sv', 'Ba': 'sv', 'Ca': 'pv', 'La': '',
    'Ti': 'pv', 'Zr': 'sv', 'Fe': 'pv', 'Mn': 'pv',
    'O': '',
}
```

### Battery Materials

```python
setups = {
    'Li': 'sv', 'Na': 'pv',
    'Fe': 'pv', 'Mn': 'pv', 'Co': '',
    'O': '', 'P': '',
}
```

### Semiconductors

```python
# Standard potentials usually sufficient
setups = {}  # Use defaults

# For high-accuracy band gaps
setups = {'Si': '', 'Ge': 'd'}
```

## Validation Checklist

Before production calculations:

1. **Check ENMAX**: Ensure ENCUT > max(ENMAX) from all POTCARs
2. **Convergence test**: Test energy convergence with ENCUT for your specific setup
3. **Compare variants**: For critical properties, compare results with different setups
4. **Check VASP wiki**: Consult official recommendations for your elements

## Running This Example

```bash
python run.py
```

The script demonstrates:
- Setting up calculations with different POTCAR variants
- Comparing energies between standard and _pv potentials
- Checking ENMAX values for proper ENCUT selection

## References

- [VASP Wiki: Available PAW potentials](https://www.vasp.at/wiki/index.php/Available_PAW_potentials)
- [VASP Wiki: Recommended potentials](https://www.vasp.at/wiki/index.php/Category:Potentials)
