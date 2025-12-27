# Configuration

## Environment Variables

| Variable | Description | Example |
|----------|-------------|---------|
| `ASE_VASP_COMMAND` | VASP executable | `mpirun -np 4 vasp_std` |
| `VASP_PP_PATH` | Pseudopotential directory | `/opt/vasp/potpaw_PBE` |
| `VASP_WORKFLOW_ENGINE` | Workflow engine | `prefect`, `dask`, `parsl` |

## Runner Configuration

### Local Runner

```python
from vasp.runners import LocalRunner

runner = LocalRunner(
    command="mpirun -np 4 vasp_std",
    pp_path="/path/to/potpaw_PBE",
    scratch_dir="/tmp/vasp",
)
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
    qos='normal',
)
```

### Kubernetes Runner

```python
from vasp import get_kubernetes_runner
KubernetesRunner = get_kubernetes_runner()

runner = KubernetesRunner(
    namespace='vasp-jobs',
    image='my-vasp-image:latest',
    cpu='4',
    memory='16Gi',
)
```

## Pseudopotential Configuration

The calculator automatically finds pseudopotentials in:

1. `VASP_PP_PATH` environment variable
2. `~/.vasp/potpaw_PBE/` default location

Directory structure:

```
$VASP_PP_PATH/
├── potpaw_PBE/
│   ├── Si/POTCAR
│   ├── Si_sv/POTCAR
│   ├── Fe/POTCAR
│   └── ...
├── potpaw_GGA/
└── potpaw_LDA/
```

## Workflow Engine

Set the workflow engine for HPC integration:

```bash
export VASP_WORKFLOW_ENGINE=prefect
```

Then recipes become workflow tasks:

```python
from vasp.recipes import relax_job

# Now decorated with @prefect.task
result = relax_job(atoms)
```
