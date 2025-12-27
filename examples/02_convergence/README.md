# 02 - Convergence Testing

Before running production calculations, you must verify that your results are converged with respect to key computational parameters.

## What You'll Learn

- Why convergence testing is essential
- How to converge ENCUT (plane-wave cutoff)
- How to converge k-point sampling
- How to choose appropriate values

## Key Parameters to Converge

### 1. ENCUT (Plane-wave Cutoff)
- Controls the size of the plane-wave basis set
- Higher ENCUT = more accurate but slower
- Rule of thumb: Use ENMAX from POTCAR × 1.3

### 2. K-points
- Controls sampling of the Brillouin zone
- Denser grids = more accurate but slower
- Metals need denser grids than insulators

## Convergence Criteria

A parameter is converged when increasing it further changes the energy by less than your required precision:
- **1 meV/atom**: Standard for most calculations
- **0.1 meV/atom**: High precision (phase stability, small energy differences)

## Running the Tests

```bash
python run.py
```

This will run a series of calculations and produce:
- `encut_convergence.png` - Energy vs ENCUT plot
- `kpoint_convergence.png` - Energy vs k-point density plot

## Typical Results for Silicon

| ENCUT (eV) | ΔE (meV/atom) |
|------------|---------------|
| 200 | 50 |
| 300 | 5 |
| 400 | 1 |
| 500 | <0.1 |

| k-grid | ΔE (meV/atom) |
|--------|---------------|
| 2×2×2 | 30 |
| 4×4×4 | 3 |
| 6×6×6 | 0.5 |
| 8×8×8 | <0.1 |

## Tips

1. Always converge ENCUT and k-points independently
2. Converge before comparing different structures
3. Similar structures may need the same parameters (error cancellation)
4. Surfaces and defects typically need higher k-points in the non-periodic direction
