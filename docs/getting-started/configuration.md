# Configuration

## Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `VASP_PP_PATH` | Pseudopotential directory | `/opt/vasp/potpaw_PBE` |
| `ASE_VASP_COMMAND` | VASP executable (optional) | `mpirun -np 4 vasp_std` |

## Pseudopotential Setup

### Obtaining POTCARs

VASP pseudopotentials (POTCARs) are **only available to licensed VASP users**. They are not open source and cannot be redistributed.

To obtain POTCARs:
1. Purchase a VASP license from [https://www.vasp.at/](https://www.vasp.at/)
2. Download the `potpaw_PBE.tgz` archive from the VASP portal
3. Extract to a permanent location

### Directory Structure

Set `VASP_PP_PATH` to point to your pseudopotential directory:

```bash
export VASP_PP_PATH=/opt/vasp/potentials
```

Expected structure:

```
$VASP_PP_PATH/
├── potpaw_PBE/           # PBE functional (most common)
│   ├── Si/
│   │   └── POTCAR
│   ├── Si_sv/            # Semi-core variant
│   │   └── POTCAR
│   ├── Fe/
│   │   └── POTCAR
│   ├── Fe_pv/            # p-valence variant
│   │   └── POTCAR
│   └── ...
├── potpaw_LDA/           # LDA functional
└── potpaw_GGA/           # Old GGA (rarely used)
```

### Selecting Variants

Use the `setups` parameter to select POTCAR variants:

```python
from vasp import Vasp

calc = Vasp(
    'my_calc',
    atoms=atoms,
    setups={
        'Fe': 'pv',   # Use Fe_pv (includes 3p electrons)
        'Li': 'sv',   # Use Li_sv (includes 1s electron)
        'O': '',      # Use standard O
    },
)
```

See [Example 19: Pseudopotential Selection](../../examples/19_pseudopotentials/) for detailed guidance on choosing variants.

## Runner Configuration

### Local Runner

```python
from vasp.runners import LocalRunner

runner = LocalRunner(
    vasp_command='vasp_std',       # VASP executable
    mpi_command='mpirun -np 4',    # MPI launcher (optional)
    background=False,              # Wait for completion
)
```

### Interactive Runner

For geometry optimizations with wavefunction reuse:

```python
from vasp.runners import InteractiveRunner

runner = InteractiveRunner(
    vasp_command='vasp_std',
    mpi_command='mpirun -np 4',
    timeout=3600,
)

# Use as context manager
with runner:
    results = runner.start('calc/', atoms)
    for step in range(100):
        if converged(results.forces):
            break
        atoms.positions = update(atoms, results.forces)
        results = runner.step(atoms)
```

### SLURM Runner

```python
from vasp.runners import SlurmRunner

runner = SlurmRunner(
    partition='compute',
    nodes=2,
    ntasks_per_node=32,
    time='4:00:00',
    account='my_project',
)
```

### Kubernetes Runner

```python
from vasp.runners import KubernetesRunner

runner = KubernetesRunner(
    namespace='vasp-jobs',
    image='my-vasp-image:latest',
    cpu='4',
    memory='16Gi',
)
```

## Quick Start

Minimal setup:

```bash
# 1. Set pseudopotential path
export VASP_PP_PATH=/path/to/potentials

# 2. Run calculation
python -c "
from ase.build import bulk
from vasp import Vasp

atoms = bulk('Si')
calc = Vasp('si_test', atoms=atoms, encut=300, kpts=(4,4,4))
energy = calc.potential_energy
print(f'Energy: {energy:.4f} eV')
"
```
