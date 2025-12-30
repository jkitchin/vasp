# VASP-ASE Interface

A modern Python interface for VASP (Vienna Ab initio Simulation Package) through ASE (Atomic Simulation Environment).

## Project Structure

```
vasp/
├── vasp/                    # Main package
│   ├── __init__.py          # Vasp calculator class
│   ├── parameters.py        # Parameter presets (VdW, DFT+U, HSE06, etc.)
│   ├── runners/             # Execution backends
│   │   ├── local.py         # LocalRunner for direct execution
│   │   ├── interactive.py   # InteractiveRunner for persistent VASP
│   │   ├── socket_io.py     # Socket server/client (i-PI protocol)
│   │   ├── slurm.py         # SLURM cluster runner
│   │   └── kubernetes.py    # Kubernetes runner
│   ├── recipes/             # Workflow recipes (quacc-style)
│   │   ├── core.py          # static_job, relax_job, double_relax_flow
│   │   ├── slabs.py         # Surface calculation recipes
│   │   ├── phonons.py       # Phonon calculation recipes
│   │   └── decorators.py    # @job, @flow, @subflow decorators
│   ├── database/            # Vector database for embeddings
│   └── tests/               # Test suite
├── docs/                    # Jupyter Book documentation
│   └── tutorials/           # 21 progressive tutorial notebooks
└── .claude/                 # Claude Code configuration
```

## Key Files

- `vasp/__init__.py` - Main `Vasp` calculator class
- `vasp/parameters.py` - Preset functions: `get_vdw_params()`, `get_ldau_params()`, `get_hybrid_params()`
- `vasp/runners/` - Execution backends: `LocalRunner`, `InteractiveRunner`, `SocketServer`, `SlurmRunner`
- `vasp/recipes/core.py` - `static_job`, `relax_job`, `double_relax_flow`
- `vasp/database/` - Vector database with per-atom embeddings

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
3. **Advanced (11-20)**: Phonons, DFT+U, HSE06, vdW, workflows, NEB, vibrations, visualization, pseudopotentials, interactive mode

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
| `/tutorial <n>` | View tutorial n (1-16) |
| `/test` | Run test suite |
| `/new-example` | Create new example |
| `/architecture` | Review codebase |
| `/lint` | Run code quality checks |
| `/build-docs` | Build Jupyter Book |
| `/status` | Project status |

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
- **vasp** - General VASP calculation help
- **job-watcher** - Job monitoring and troubleshooting
- **troubleshoot** - Error diagnosis and fixes
