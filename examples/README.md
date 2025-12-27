# VASP Examples

A progressive collection of examples demonstrating the vasp-ase interface, from basic calculations to advanced workflows.

## Prerequisites

- VASP installed and accessible (set `ASE_VASP_COMMAND` or configure runner)
- ASE and this package installed
- VASP pseudopotentials (set `VASP_PP_PATH`)

## Directory Structure

### Beginner (Start Here)

| Directory | Topic | Time | Description |
|-----------|-------|------|-------------|
| `01_getting_started/` | First calculation | 1 min | Energy of a silicon crystal |
| `02_convergence/` | Convergence testing | 5 min | ENCUT and k-point convergence |
| `03_relaxation/` | Structure optimization | 2 min | Relaxing atomic positions and cell |
| `04_equation_of_state/` | Bulk modulus | 3 min | Fitting energy vs volume |
| `05_density_of_states/` | Electronic structure | 2 min | Total and projected DOS |

### Intermediate

| Directory | Topic | Time | Description |
|-----------|-------|------|-------------|
| `06_band_structure/` | Band structure | 3 min | Band structure along high-symmetry path |
| `07_magnetism/` | Spin polarization | 3 min | Magnetic moments in Fe and NiO |
| `08_surfaces/` | Slab calculations | 5 min | Surface energy of Cu(111) |
| `09_adsorption/` | Molecules on surfaces | 5 min | CO on Pt(111) |
| `10_reactions/` | Reaction energetics | 5 min | Dissociation energy and barriers |

### Advanced

| Directory | Topic | Time | Description |
|-----------|-------|------|-------------|
| `11_phonons/` | Vibrational properties | 10 min | Phonon dispersion with Phonopy |
| `12_dft_plus_u/` | Strongly correlated | 3 min | DFT+U for FeO |
| `13_hybrid_functionals/` | HSE06 | 10 min | Band gap with hybrid functional |
| `14_van_der_waals/` | Dispersion corrections | 3 min | Graphite interlayer binding |
| `15_workflows/` | Automated pipelines | 5 min | Using recipes for high-throughput |
| `16_neb/` | Nudged elastic band | 15 min | Transition states and reaction barriers |

## Running Examples

Each directory contains:
- `run.py` - Main executable script
- `README.md` - Detailed explanation
- Optional: `analyze.py` - Post-processing script

To run an example:
```bash
cd examples/01_getting_started
python run.py
```

## Tips for Fast Calculations

These examples use small systems and coarse settings for speed:
- Low ENCUT (300-400 eV)
- Sparse k-point grids
- Minimal supercells

For production calculations, increase these parameters based on convergence testing.

## Environment Setup

```bash
# Set VASP command (adjust for your system)
export ASE_VASP_COMMAND="mpirun -np 4 vasp_std"

# Set pseudopotential path
export VASP_PP_PATH=/path/to/potpaw_PBE

# Or use in Python:
from vasp import Vasp
from vasp.runners import LocalRunner

runner = LocalRunner(
    command="mpirun -np 4 vasp_std",
    pp_path="/path/to/potpaw_PBE"
)
```
