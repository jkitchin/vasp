# VASP-ASE Interface

[![Tests](https://github.com/jkitchin/vasp/actions/workflows/tests.yml/badge.svg)](https://github.com/jkitchin/vasp/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/jkitchin/vasp/branch/master/graph/badge.svg)](https://codecov.io/gh/jkitchin/vasp)
[![Documentation](https://img.shields.io/badge/docs-online-blue)](https://kitchingroup.cheme.cmu.edu/vasp/)

> **⚠️ ALPHA SOFTWARE**: This project is under active development and not ready for production use. APIs may change without notice. This notice will be removed when the project reaches stable release.

A modern Python interface for [VASP](https://www.vasp.at/) (Vienna Ab initio Simulation Package) through [ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment).

## Features

- **Modern Calculator**: Clean ASE-compatible interface with type hints
- **Parameter Presets**: Built-in support for DFT+U, hybrid functionals, van der Waals, spin-orbit coupling
- **Multiple Runners**: Local, SLURM, Kubernetes execution backends
- **Workflow Recipes**: quacc-style `@job` and `@flow` decorators for complex workflows
- **Analysis Tools**: DOS, band structure, Bader charges, and more

## Installation

```bash
pip install vasp-ase
```

For development:
```bash
git clone https://github.com/jkitchin/vasp.git
cd vasp
pip install -e ".[dev]"
```

## Quick Start

```python
from ase.build import bulk
from vasp import Vasp

# Create atoms
atoms = bulk('Cu', 'fcc', a=3.6)

# Run calculation
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
    directory='cu-bulk'
)

energy = atoms.get_potential_energy()
print(f"Energy: {energy:.3f} eV")
```

## Parameter Presets

```python
from vasp.parameters import get_vdw_params, get_ldau_params, get_hybrid_params

# Van der Waals (D3-BJ)
vdw = get_vdw_params('d3bj')

# DFT+U for transition metals
ldau = get_ldau_params(['Fe', 'O'], {'Fe': {'U': 4.0, 'J': 0.0}})

# Hybrid functionals
hse = get_hybrid_params('hse06')
```

## Workflow Recipes

```python
from vasp.recipes.core import static_job, relax_job, double_relax_flow

# Single relaxation
result = relax_job(atoms, kpts=(4, 4, 4))

# Double relaxation workflow
result = double_relax_flow(atoms)
```

## Runners

```python
from vasp.runners import LocalRunner, SlurmRunner

# Local execution
runner = LocalRunner(vasp_command='mpirun -np 4 vasp_std')

# SLURM cluster
runner = SlurmRunner(
    partition='compute',
    nodes=2,
    ntasks_per_node=24,
    time='2:00:00'
)
```

## Tutorials

The `docs/tutorials/` directory contains 21 progressive tutorials as Jupyter notebooks:

1. **Beginner (01-05)**: Energy, convergence, relaxation, EOS, DOS
2. **Intermediate (06-10)**: Bands, magnetism, surfaces, adsorption, reactions
3. **Advanced (11-21)**: Phonons, DFT+U, HSE06, vdW, workflows, NEB, vibrations, visualization, pseudopotentials, interactive mode, cluster expansion

Run tutorials with:
```bash
make tutorial T=01    # Run a single tutorial
make tutorials        # Run all tutorials
```

## Claude Code Integration

Install global Claude Code skills for VASP assistance:

```bash
pip install vasp-ase
vasp-claude install
```

Then use commands like `/vasp-help`, `/vasp-watch-job`, `/vasp-fix-job` from any project.

## Documentation

Full documentation available at https://kitchingroup.cheme.cmu.edu/vasp/

## Related Projects

This project is part of a growing ecosystem of LLM agents for atomistic simulation:

| Project | Description | Reference |
|---------|-------------|-----------|
| [DREAMS](https://arxiv.org/abs/2507.14267) | Multi-agent DFT simulation framework | Wang & Viswanathan, 2025 |
| [AtomAgents](https://github.com/lamm-mit/AtomAgents) | Physics-aware alloy design with LAMMPS | Ghafarollahi & Buehler, 2025 |
| [LangSim](https://github.com/jan-janssen/LangSim) | LLM interface for atomistic simulation with pyiron | Janssen, 2024 |
| [El Agente](https://doi.org/10.1016/j.matt.2025.102263) | Autonomous quantum chemistry workflows | Zou & Aspuru-Guzik, 2025 |

See [references.bib](references.bib) for full citations.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please read the contributing guidelines and submit pull requests.

## Citation

If you use this package, please cite:
- VASP: https://www.vasp.at/
- ASE: https://wiki.fysik.dtu.dk/ase/
