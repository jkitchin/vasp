# 03 - Structure Relaxation

Most structures from experiment or initial guesses are not at their theoretical equilibrium. Geometry optimization (relaxation) finds the minimum energy configuration.

## What You'll Learn

- How to relax atomic positions
- How to relax the unit cell (volume and shape)
- Understanding ISIF parameter
- Convergence criteria for forces

## Types of Relaxation

| ISIF | Ions | Cell Shape | Cell Volume | Use Case |
|------|------|------------|-------------|----------|
| 0 | No | No | No | Single-point energy |
| 2 | Yes | No | No | Surfaces, defects |
| 3 | Yes | Yes | Yes | Bulk equilibrium |
| 4 | Yes | Yes | No | Fixed volume relaxation |

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `ibrion` | 2 | Conjugate gradient optimizer |
| `isif` | 2 or 3 | What to relax (see above) |
| `nsw` | 100 | Maximum ionic steps |
| `ediffg` | -0.02 | Force convergence (negative = force criterion in eV/Å) |

## The Calculations

This example performs two relaxations:

1. **Ionic relaxation (ISIF=2)**: Optimize atomic positions at fixed cell
2. **Full relaxation (ISIF=3)**: Optimize both atoms and cell parameters

## Expected Results

Starting from a slightly distorted silicon structure:
- Forces reduce from ~0.5 eV/Å to <0.02 eV/Å
- Lattice constant relaxes to ~5.47 Å (PBE overestimates by ~1%)

## Tips

1. Always start from a reasonable initial structure
2. Use ISIF=2 for surfaces (don't relax the vacuum)
3. For metals, use smaller EDIFFG or more steps
4. Check that forces are actually converged in OUTCAR
