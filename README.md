# VASP-ASE Interface

A modern Python interface for VASP (Vienna Ab initio Simulation Package) through ASE (Atomic Simulation Environment).

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

The `examples/` directory contains 16 progressive tutorials:

1. **Beginner (01-05)**: Energy, convergence, relaxation, EOS, DOS
2. **Intermediate (06-10)**: Bands, magnetism, surfaces, adsorption, reactions
3. **Advanced (11-16)**: Phonons, DFT+U, HSE06, vdW, workflows, NEB

## Claude Code Integration

Install global Claude Code skills for VASP assistance:

```bash
pip install vasp-ase
vasp-claude install
```

Then use commands like `/vasp-help`, `/vasp-watch-job`, `/vasp-fix-job` from any project.

## Documentation

Full documentation available at the [project docs](https://github.com/jkitchin/vasp).

## License

MIT License - see [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome! Please read the contributing guidelines and submit pull requests.

## Citation

If you use this package, please cite:
- VASP: https://www.vasp.at/
- ASE: https://wiki.fysik.dtu.dk/ase/
