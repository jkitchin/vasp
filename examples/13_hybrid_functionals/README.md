# 13 - Hybrid Functionals

Calculate accurate band gaps using HSE06 hybrid functional.

## Key Concepts

- **Exact Exchange**: Part of HF exchange mixed with DFT
- **HSE06**: Screened hybrid, most popular for solids
- **Band Gap**: Major improvement over GGA
- **Computational Cost**: 10-100x more than GGA

## Requirements

Hybrid calculations require significant computational resources!

## VASP Parameters

### HSE06 Setup

```python
from vasp.parameters import get_hybrid_params

hse = get_hybrid_params('hse06')

calc = Vasp(
    atoms=atoms,
    encut=400,
    kpts=(4, 4, 4),  # Fewer than GGA
    ismear=0,
    sigma=0.05,
    **hse,
)
```

### Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `LHFCALC` | True | Enable hybrid |
| `HFSCREEN` | 0.2 | Screening (HSE) |
| `AEXX` | 0.25 | Fraction of HF |
| `ALGO` | Damped | Mixer for hybrids |

## Available Hybrids

| Functional | AEXX | Screening | Use Case |
|------------|------|-----------|----------|
| PBE0 | 0.25 | None | Small molecules |
| HSE06 | 0.25 | 0.2 Å⁻¹ | Solids (standard) |
| HSE03 | 0.25 | 0.3 Å⁻¹ | Solids (less screening) |
| B3LYP | 0.20 | None | Molecules |

## Physical Background

### Why Hybrids Work

Standard DFT:
- Self-interaction error overdelocalizes electrons
- Band gaps underestimated (often by 30-50%)

Hybrid functionals:
- Mix exact exchange: E_xc = a×E_x^HF + (1-a)×E_x^DFT + E_c^DFT
- Reduce self-interaction error
- Improve band gap and localization

### Band Gap Improvement

| Material | GGA | HSE06 | Exp |
|----------|-----|-------|-----|
| Si | 0.6 eV | 1.1 eV | 1.12 eV |
| GaAs | 0.4 eV | 1.4 eV | 1.42 eV |
| MgO | 4.5 eV | 7.1 eV | 7.8 eV |
| ZnO | 0.7 eV | 2.5 eV | 3.3 eV |

## Output

```
Silicon with HSE06:
  Total energy: -10.87 eV
  Estimated band gap: ~1.1 eV
  Experimental: 1.12 eV
```

## Computational Tips

1. **Start from PBE**:
   ```python
   # First: PBE relaxation
   # Then: HSE static with ISTART=1
   ```

2. **Reduce cost**:
   ```
   PRECFOCK = Fast    # Coarser HF grid
   NKRED = 2          # Halve HF k-mesh
   ```

3. **K-points**: Hybrids need fewer k-points than GGA

4. **Memory**: ALGO = All if memory limited

## When to Use

**Use HSE06:**
- Accurate band gaps needed
- Defect level alignment
- Charge localization
- Band offsets at interfaces

**Use GGA:**
- Geometry optimization
- Initial screening
- Metals
- Large systems

## Applications

- Semiconductor band gaps
- Defect physics
- Interface band alignment
- Photovoltaic materials

## Next Steps

Try `14_van_der_waals/` for dispersion-corrected calculations.
