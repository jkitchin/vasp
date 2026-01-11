# VASP-ASE Interface

A modern Python interface for VASP (Vienna Ab initio Simulation Package) through ASE (Atomic Simulation Environment).

**Package:** vasp-ase v2.0.0
**Python:** ≥3.9 (supports 3.9-3.13)
**License:** MIT

## Project Structure

```
vasp/
├── vasp/                          # Main Python package
│   ├── __init__.py               # Package exports and main Vasp class
│   ├── calculator.py             # Core Vasp calculator implementation
│   ├── parameters.py             # Parameter presets and validation
│   ├── exceptions.py             # Exception hierarchy (async patterns)
│   ├── cli.py                    # CLI tools and commands
│   ├── mixins/                   # Modular functionality mixins
│   │   ├── io.py                 # File I/O (INCAR, POSCAR, KPOINTS, POTCAR)
│   │   ├── electronic.py         # Electronic structure (bands, DOS)
│   │   ├── analysis.py           # Post-processing (Bader, unfolding)
│   │   └── dynamics.py           # NEB and vibrational calculations
│   ├── runners/                  # Execution backends
│   │   ├── base.py               # Abstract Runner interface
│   │   ├── local.py              # LocalRunner for direct execution
│   │   ├── mock.py               # MockRunner for testing without VASP
│   │   ├── interactive.py        # InteractiveRunner for persistent VASP
│   │   ├── socket_io.py          # Socket server/client (i-PI protocol)
│   │   ├── slurm.py              # SLURM cluster runner
│   │   └── kubernetes.py         # Kubernetes runner
│   ├── recipes/                  # Workflow recipes (quacc-style)
│   │   ├── core.py               # static_job, relax_job, double_relax_flow
│   │   ├── slabs.py              # Surface calculation recipes
│   │   ├── phonons.py            # Phonon calculation recipes
│   │   └── decorators.py         # @job, @flow, @subflow decorators
│   ├── database/                 # Vector database for embeddings
│   │   ├── vector_db.py          # VectorAtomDatabase with libSQL
│   │   └── embedders.py          # Atomic environment embedders
│   └── tests/                    # Test suite
│       ├── conftest.py           # Pytest fixtures
│       ├── fixtures/             # Mock VASP output data
│       └── test_*.py             # 11 test modules
├── docs/                         # Jupyter Book documentation
│   ├── _config.yml               # Jupyter Book configuration
│   ├── _toc.yml                  # Table of contents
│   ├── api/                      # API documentation
│   ├── guide/                    # User guides
│   └── tutorials/                # 22 interactive notebooks
├── vscode-vasp-syntax/           # VS Code syntax highlighting extension
├── .claude/                      # Claude Code configuration
│   ├── settings.json             # SessionStart hook
│   ├── commands/                 # 13 slash commands
│   └── skills/                   # 3 AI assistant skills
├── pyproject.toml                # Package configuration
├── Makefile                      # Development targets
└── *.md                          # Project documentation
```

## Key Modules

| Module | Purpose | Key Exports |
|--------|---------|-------------|
| `vasp/__init__.py` | Main package entry | `Vasp` |
| `vasp/calculator.py` | Calculator implementation | `Vasp`, `CalculationResult` |
| `vasp/parameters.py` | Parameter presets | `get_vdw_params()`, `get_ldau_params()`, `get_hybrid_params()`, `get_soc_params()`, `get_md_params()`, `get_phonon_params()`, `HubbardU` |
| `vasp/exceptions.py` | Async-friendly exceptions | `VaspSubmitted`, `VaspQueued`, `VaspRunning`, `VaspNotFinished`, `VaspNotConverged`, `VaspError` |
| `vasp/cli.py` | CLI tools | `main()`, `vaspsum()`, `vasp_debug()` |
| `vasp/runners/` | Execution backends | `LocalRunner`, `MockRunner`, `SlurmRunner`, `KubernetesRunner`, `InteractiveRunner`, `SocketServer` |
| `vasp/recipes/core.py` | Core workflows | `static_job()`, `relax_job()`, `double_relax_flow()` |
| `vasp/recipes/phonons.py` | Phonon workflows | `phonon_flow()`, `dfpt_phonon_job()` |
| `vasp/recipes/slabs.py` | Surface workflows | `slab_static_job()`, `slab_relax_job()`, `bulk_to_slabs_flow()` |
| `vasp/database/vector_db.py` | Structure database | `VectorAtomDatabase`, `StructureRecord` |

## Development Commands

### Makefile Targets

```bash
make install          # Install in editable mode
make install-dev      # Install with dev dependencies
make test             # Run pytest suite
make test-coverage    # Generate coverage report
make lint             # Run ruff linter
make typecheck        # Run mypy type checker
make format           # Auto-format code
make check            # Run all checks (lint + typecheck + test)
make docs-build       # Build Jupyter Book
make docs-serve       # Build and serve docs locally
make tutorials        # Run all tutorial notebooks
make tutorial T=<n>   # Run single tutorial
```

### Direct Commands

```bash
pytest                        # Run tests
pytest --cov=vasp             # Run tests with coverage
ruff check .                  # Lint code
ruff format .                 # Format code
mypy vasp                     # Type check
jupyter-book build docs       # Build documentation
```

## Runners

| Runner | Purpose | Use Case |
|--------|---------|----------|
| `LocalRunner` | Direct local execution | Development, small jobs |
| `MockRunner` | Simulated VASP | Unit tests, workflow design |
| `SlurmRunner` | SLURM cluster submission | HPC clusters |
| `KubernetesRunner` | Kubernetes orchestration | Cloud computing |
| `InteractiveRunner` | Persistent VASP session | Low-latency repeated calcs |
| `SocketServer/Client` | i-PI socket protocol | MD coupling |

All runners implement the abstract `Runner` interface from `runners/base.py`.

## Recipes

### Core (`recipes/core.py`)
- `static_job()` - Single-point energy
- `relax_job()` - Ionic/cell relaxation
- `double_relax_flow()` - Sequential ionic then cell relaxation
- `static_from_relax_job()` - Static after relaxation
- `relax_and_static_flow()` - Combined workflow

### Phonons (`recipes/phonons.py`)
- `phonon_flow()` - Complete phonon calculation
- `phonon_displacement_job()` - Single displacement
- `dfpt_phonon_job()` - DFPT phonons

### Surfaces (`recipes/slabs.py`)
- `slab_static_job()` - Static on slab
- `slab_relax_job()` - Relax slab
- `bulk_to_slabs_flow()` - Complete slab workflow
- `calculate_surface_energy()` - Surface energy

### Decorators (`recipes/decorators.py`)
- `@job` - Single calculation job
- `@flow` - Workflow of jobs
- `@subflow` - Nested workflow

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

### Parameter Presets

```python
# DFT+U
from vasp.parameters import get_ldau_params, HubbardU
params = get_ldau_params(['Fe', 'O'], {'Fe': HubbardU(u=4.0)})

# Hybrid functionals
from vasp.parameters import get_hybrid_params
params = get_hybrid_params('hse06')  # Also: 'pbe0', 'b3lyp'

# Van der Waals
from vasp.parameters import get_vdw_params
params = get_vdw_params('d3bj')  # Also: 'd3', 'd2', 'optb88', 'ts'

# Spin-orbit coupling
from vasp.parameters import get_soc_params
params = get_soc_params()

# Molecular dynamics
from vasp.parameters import get_md_params
params = get_md_params(ensemble='nvt', temperature=300)

# Phonons
from vasp.parameters import get_phonon_params
params = get_phonon_params()
```

## Tutorial Notebooks (22 total)

### Beginner (00-05)
- **00** System info and environment verification
- **01** Getting started - first calculation
- **02** Convergence - ENCUT and k-points
- **03** Relaxation - structure optimization
- **04** Equation of state - bulk modulus
- **05** Density of states

### Intermediate (06-10)
- **06** Band structure
- **07** Magnetism - spin-polarized DFT
- **08** Surfaces - slab calculations
- **09** Adsorption - molecules on surfaces
- **10** Reactions - energetics

### Advanced (11-21)
- **11** Phonons - dispersions and DOS
- **12** DFT+U - correlated systems
- **13** Hybrid functionals - HSE06, PBE0
- **14** Van der Waals - dispersion corrections
- **15** Workflows - automated pipelines
- **16** NEB - transition states
- **17** Vibrations - IR/Raman
- **18** 3D visualization
- **19** Pseudopotentials - POTCAR setup
- **20** Interactive mode - persistent VASP
- **21** Cluster expansion - alloy thermodynamics

## Testing

Tests use `MockRunner` to run without VASP:

```python
from vasp.runners import MockRunner, MockResults
import numpy as np

mock = MockResults(energy=-10.5, forces=np.zeros((2, 3)))
runner = MockRunner(results=mock)
calc = Vasp(atoms=atoms, runner=runner, xc='PBE')
```

### Test Modules
- `test_calculator.py` - Calculator initialization and execution
- `test_runners.py` - All runner implementations
- `test_parameters.py` - Parameter preset validation
- `test_recipes.py` - Workflow job and flow execution
- `test_database.py` - Database operations
- `test_vector_database.py` - Vector similarity search
- `test_analysis.py` - Post-processing methods
- `test_mixins.py` - Mixin functionality
- `test_socket_io.py` - Socket communication
- `test_interactive.py` - Interactive runner
- `test_exceptions.py` - Exception handling

## Async/Exception-Based Patterns

The package uses exceptions for async job status signaling:

```python
from vasp.exceptions import VaspSubmitted, VaspQueued, VaspRunning, VaspNotFinished

try:
    result = calc.get_potential_energy()
except VaspSubmitted as e:
    print(f"Job submitted: {e.job_id}")
except VaspQueued as e:
    print(f"Job queued, position: {e.position}")
except VaspRunning as e:
    print(f"Job running: {e.progress}%")
except VaspNotFinished:
    print("Job not yet complete")
```

Alternative: Use `CalculationResult` for exception-free status checking.

## Code Style

- Use type hints for all function signatures
- Follow existing patterns in the codebase
- Add tests for new features (use `MockRunner`)
- Update documentation for user-facing changes
- Run `make check` before committing

## Claude Code Integration

### Global Installation

Install Claude Code skills globally (works from any project):

```bash
vasp-claude install    # Install skills
vasp-claude status     # Check installation
vasp-claude uninstall  # Remove skills
```

### Project Commands

Use these slash commands when working in this repository:

| Command | Description |
|---------|-------------|
| `/docs` | Open documentation |
| `/examples` | List all examples |
| `/tutorial <n>` | View tutorial n (0-21) |
| `/notebooks` | Open and explore notebooks |
| `/test` | Run test suite |
| `/lint` | Run code quality checks |
| `/build-docs` | Build Jupyter Book |
| `/status` | Project status |
| `/architecture` | Review codebase |
| `/new-example` | Create new example |

### Job Monitoring Commands

| Command | Description |
|---------|-------------|
| `/watch-job <dir>` | Monitor VASP job status |
| `/fix-job <dir>` | Diagnose and fix failed job |
| `/vasp-help <topic>` | Parameter reference |

### Global Commands (after `vasp-claude install`)

These work from any project:

| Command | Description |
|---------|-------------|
| `/vasp-help <topic>` | VASP parameter help |
| `/vasp-watch-job <dir>` | Monitor running job |
| `/vasp-fix-job <dir>` | Auto-fix failed job |
| `/vasp-examples` | List tutorials |
| `/vasp-tutorial <n>` | View tutorial |

### Skills

Claude automatically uses these skills:
- **vasp-calculation** - General VASP calculation help
- **job-watcher** - Job monitoring and troubleshooting
- **troubleshoot** - Error diagnosis and fixes

## Dependencies

### Required
- numpy ≥1.22
- ase ≥3.22
- scipy ≥1.8

### Optional Groups
- `dev` - pytest, pytest-cov, mypy, ruff, build, twine
- `kubernetes` - kubernetes ≥28.0
- `phonons` - phonopy ≥2.20
- `plotting` - matplotlib ≥3.5
- `docs` - jupyter-book, sphinx, myst-nb, sphinxcontrib-bibtex

Install optional dependencies:
```bash
pip install -e ".[dev,phonons,plotting]"
```

## VS Code Extension

The `vscode-vasp-syntax/` directory contains a VS Code extension for VASP file syntax highlighting:
- INCAR files
- POSCAR/CONTCAR files
- KPOINTS files
- OUTCAR files

Build and install with `npm install && npm run compile` from the extension directory.
