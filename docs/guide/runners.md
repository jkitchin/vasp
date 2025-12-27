# Runners

Runners handle the execution of VASP calculations on different backends.

## Local Runner

Runs VASP on the local machine:

```python
from vasp.runners import LocalRunner

runner = LocalRunner(
    command="mpirun -np 4 vasp_std",
    pp_path="/path/to/potpaw_PBE",
    scratch_dir="/tmp/vasp",
)

calc = Vasp(atoms=atoms, runner=runner, ...)
```

## SLURM Runner

Submits jobs to SLURM scheduler:

```python
from vasp import get_slurm_runner
SlurmRunner = get_slurm_runner()

runner = SlurmRunner(
    partition='compute',
    nodes=2,
    ntasks_per_node=32,
    time='4:00:00',
    account='my_project',
)
```

### Non-blocking Execution

```python
from vasp import VaspSubmitted, VaspQueued, VaspRunning

try:
    energy = calc.potential_energy
except VaspSubmitted as e:
    print(f"Job submitted: {e.jobid}")
except VaspQueued as e:
    print(f"Job queued: {e.jobid}")
except VaspRunning as e:
    print(f"Job running: {e.jobid}")
```

## Kubernetes Runner

Runs VASP in Kubernetes pods:

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

## Mock Runner

For testing without VASP:

```python
from vasp.runners import MockRunner, MockResults

mock = MockResults(energy=-10.5)
runner = MockRunner(results=mock)

calc = Vasp(atoms=atoms, runner=runner, ...)
energy = calc.potential_energy  # Returns -10.5
```

## Custom Runners

Implement the `Runner` base class:

```python
from vasp.runners import Runner, JobStatus

class MyRunner(Runner):
    def submit(self, directory):
        # Submit job
        return job_id

    def status(self, job_id):
        # Check status
        return JobStatus(state=JobState.RUNNING)
```
