# VASP-ASE Interface

A modern Python interface for VASP (Vienna Ab initio Simulation Package) through ASE (Atomic Simulation Environment).

## Features

- **Full VASP Parameter Support**: All VASP parameters via keyword arguments
- **Automatic Input/Output**: Generate inputs, parse outputs automatically
- **Pluggable Runners**: Local, SLURM, Kubernetes execution backends
- **Non-blocking Execution**: Exception-based job status signaling
- **Workflow Integration**: Prefect, Dask, Parsl, Covalent support
- **Parameter Presets**: Common setups for VdW, DFT+U, HSE06, etc.
- **Recipe System**: Reusable calculation patterns (quacc-style)

## Quick Example

```python
from ase.build import bulk
from vasp import Vasp

# Create silicon structure
atoms = bulk('Si', 'diamond', a=5.43)

# Run VASP calculation
calc = Vasp(
    atoms=atoms,
    xc='PBE',
    encut=400,
    kpts=(8, 8, 8),
)

energy = calc.potential_energy
print(f"Energy: {energy:.3f} eV")
```

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

## Documentation

```{tableofcontents}
```

## License

GPL-3.0-or-later
