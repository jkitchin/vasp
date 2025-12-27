# 01 - Getting Started: Your First VASP Calculation

This example demonstrates the most basic VASP calculation: computing the total energy of bulk silicon.

## What You'll Learn

- How to set up a basic VASP calculation
- Creating atoms with ASE
- Understanding key VASP parameters
- Reading results

## The Calculation

We calculate the total energy of silicon in the diamond structure using:
- PBE exchange-correlation functional
- Plane-wave cutoff of 300 eV
- 4×4×4 k-point grid

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `xc` | PBE | Exchange-correlation functional |
| `encut` | 300 | Plane-wave cutoff in eV |
| `kpts` | (4,4,4) | k-point mesh |
| `ismear` | 1 | Methfessel-Paxton smearing (for metals/semiconductors) |
| `sigma` | 0.1 | Smearing width in eV |

## Expected Output

- Total energy: approximately -10.8 eV (2 atoms)
- Runtime: ~30-60 seconds on 4 cores

## Files Generated

After running, you'll find:
- `INCAR` - Input parameters
- `POSCAR` - Crystal structure
- `POTCAR` - Pseudopotentials (linked)
- `KPOINTS` - k-point mesh
- `OUTCAR` - Detailed output
- `vasprun.xml` - XML output for parsing

## Next Steps

Once this runs successfully, proceed to `02_convergence/` to learn about convergence testing.
