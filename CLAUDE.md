# VASP-ASE Interface

A modern Python interface for VASP (Vienna Ab initio Simulation Package) through ASE (Atomic Simulation Environment).

## Project Structure

```
vasp/
├── vasp/                    # Main package
│   ├── __init__.py          # Vasp calculator class
│   ├── parameters.py        # Parameter presets (VdW, DFT+U, HSE06, etc.)
│   ├── runners.py           # Execution backends (Local, SLURM, Kubernetes)
│   ├── recipes/             # Workflow recipes (quacc-style)
│   │   ├── core.py          # static_job, relax_job, double_relax_flow
│   │   ├── slabs.py         # Surface calculation recipes
│   │   ├── phonons.py       # Phonon calculation recipes
│   │   └── decorators.py    # @job, @flow, @subflow decorators
│   └── tests/               # Test suite
├── examples/                # 16 progressive tutorials
│   ├── 01_getting_started/  # Basic energy calculation
│   ├── ...
│   └── 16_neb/              # Nudged elastic band
├── docs/                    # Jupyter Book documentation
└── .claude/                 # Claude Code configuration
```

## Key Files

- `vasp/__init__.py` - Main `Vasp` calculator class
- `vasp/parameters.py` - Preset functions: `get_vdw_params()`, `get_ldau_params()`, `get_hybrid_params()`
- `vasp/runners.py` - `LocalRunner`, `MockRunner`, `SlurmRunner`, `KubernetesRunner`
- `vasp/recipes/core.py` - `static_job`, `relax_job`, `double_relax_flow`

## Development Commands

```bash
# Run tests
pytest

# Run tests with coverage
pytest --cov=vasp

# Build documentation
jupyter-book build docs

# Lint code
ruff check .
```

## VASP Parameter Quick Reference

### Common Parameters
- `xc` - Exchange-correlation: 'PBE', 'LDA', 'PW91'
- `encut` - Plane-wave cutoff (eV)
- `kpts` - K-point mesh: (8, 8, 8) or path for bands
- `ismear` - Smearing: -5 (tetrahedron), 0 (Gaussian), 1 (MP)
- `sigma` - Smearing width (eV)

### Relaxation
- `ibrion` - Optimizer: 1 (quasi-Newton), 2 (CG)
- `isif` - What to relax: 2 (ions), 3 (ions+cell)
- `nsw` - Max ionic steps
- `ediffg` - Force convergence (negative = forces in eV/Å)

### Spin/Magnetism
- `ispin` - 1 (non-spin), 2 (spin-polarized)
- `magmom` - Initial magnetic moments

### DFT+U
- Use `get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})`

### Hybrid Functionals
- Use `get_hybrid_params('hse06')`

### Van der Waals
- Use `get_vdw_params('d3bj')` for D3-BJ correction

## Tutorial Progression

1. **Beginner (01-05)**: Energy, convergence, relaxation, EOS, DOS
2. **Intermediate (06-10)**: Bands, magnetism, surfaces, adsorption, reactions
3. **Advanced (11-16)**: Phonons, DFT+U, HSE06, vdW, workflows, NEB

## Code Style

- Use type hints
- Follow existing patterns in the codebase
- Add tests for new features
- Update documentation for user-facing changes

## Testing

Tests use `MockRunner` to run without VASP:

```python
from vasp.runners import MockRunner, MockResults

mock = MockResults(energy=-10.5, forces=np.zeros((2, 3)))
runner = MockRunner(results=mock)
calc = Vasp(atoms=atoms, runner=runner, ...)
```

## Claude Code Commands

Use these slash commands for common tasks:
- `/docs` - Open documentation
- `/examples` - List available examples
- `/tutorial <n>` - View a specific tutorial
- `/test` - Run the test suite
- `/new-example` - Create a new example
