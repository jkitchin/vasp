# Installation

## Requirements

- Python 3.10 or later
- VASP installed and accessible
- VASP pseudopotentials

## Install from PyPI

```bash
pip install vasp-ase
```

## Install from Source

```bash
git clone https://github.com/jkitchin/vasp.git
cd vasp
pip install -e .
```

## Install with Optional Dependencies

```bash
# Development tools
pip install -e ".[dev]"

# Phonon calculations
pip install -e ".[phonons]"

# Kubernetes runner
pip install -e ".[kubernetes]"

# All optional dependencies
pip install -e ".[all]"
```

## Configure VASP

Set environment variables:

```bash
# VASP executable
export ASE_VASP_COMMAND="mpirun -np 4 vasp_std"

# Pseudopotential path
export VASP_PP_PATH=/path/to/potpaw_PBE
```

Or configure in Python:

```python
from vasp import Vasp
from vasp.runners import LocalRunner

runner = LocalRunner(
    command="mpirun -np 4 vasp_std",
    pp_path="/path/to/potpaw_PBE"
)

calc = Vasp(atoms=atoms, runner=runner, ...)
```

## Verify Installation

```python
from vasp import Vasp, __version__
print(f"vasp-ase version: {__version__}")
```
