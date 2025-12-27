# 12 - DFT+U

Apply Hubbard U corrections for strongly correlated systems using FeO as an example.

## Key Concepts

- **Hubbard U**: On-site Coulomb interaction for localized electrons
- **Self-Interaction Error**: DFT's overdelocalization of d/f electrons
- **Mott Insulators**: Materials that DFT incorrectly predicts as metals
- **Dudarev Formulation**: Simplified U_eff = U - J approach

## VASP Parameters

### DFT+U Setup

```python
from vasp.parameters import get_ldau_params, HubbardU

ldau = get_ldau_params(
    symbols=['Fe', 'O'],
    u_values={'Fe': HubbardU(u=4.0, j=0.0)},
)

calc = Vasp(
    atoms=feo,
    xc='PBE',
    encut=500,
    kpts=(8, 8, 8),
    ispin=2,
    magmom=[4.0, 0.0],
    lorbit=11,
    **ldau,
)
```

### Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `LDAU` | .TRUE. | Enable DFT+U |
| `LDAUTYPE` | 2 | Dudarev formulation |
| `LDAUL` | [2, -1] | d-orbitals for Fe, none for O |
| `LDAUU` | [4.0, 0.0] | U values |
| `LDAUJ` | [0.0, 0.0] | J values |
| `LMAXMIX` | 4 | For d-electrons |

## Physical Background

### Why DFT Fails

Standard DFT (LDA/GGA):
- Overdelocalizes d/f electrons
- Underestimates band gaps
- Incorrectly predicts metals for Mott insulators

### Hubbard Model

The DFT+U energy correction:

**Dudarev (LDAUTYPE=2):**
```
E_U = (U-J)/2 × Σ_m [n_m - n_m²]
```

This penalizes fractional occupation, favoring integer occupations.

### Effect of U

| Property | DFT | DFT+U |
|----------|-----|-------|
| Localization | Poor | Good |
| Band gap | Underestimated | Improved |
| Magnetic moment | Underestimated | Improved |
| Metal prediction | Sometimes wrong | Correct |

## Output

```
Comparison:
  Method          Fe moment (μB)     Behavior
  Standard DFT    3.12               Metallic
  DFT+U (U=4)     3.58               Insulating
  Experimental    ~3.6               Insulating
```

## Typical U Values

| Element | System | U (eV) |
|---------|--------|--------|
| Fe | Oxides | 4-5 |
| Co | Oxides | 3-4 |
| Ni | Oxides | 5-7 |
| Mn | Oxides | 3-4 |
| Ce | f-electrons | 5-6 |
| U | Actinides | 4-5 |

## How to Determine U

1. **Linear Response**: Self-consistent U from constrained DFT
2. **Empirical Fitting**: Match experimental properties
3. **Literature Values**: Use established values

## Tips

1. **LMAXMIX**: Set to 4 for d-electrons, 6 for f-electrons
2. **Converge U**: Test U values for your property of interest
3. **Consistency**: Use same U for all calculations in a study
4. **Starting Occupation**: May need to adjust MAGMOM

## Common Issues

- **Multiple Minima**: Different starting configurations may give different solutions
- **U Dependence**: Some properties are very U-sensitive
- **Transferability**: U values may not transfer between systems

## Applications

- Transition metal oxides
- Lanthanide/actinide compounds
- Magnetic insulators
- Battery materials (LiFePO4, etc.)

## Next Steps

Try `13_hybrid_functionals/` for HSE06 calculations.
