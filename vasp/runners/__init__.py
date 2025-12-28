"""Runners for executing VASP calculations.

Runners handle the execution of VASP jobs, whether locally,
on a cluster scheduler (SLURM), or on Kubernetes. They provide
a uniform interface for submitting jobs and checking status.

Example:
    >>> from vasp.runners import LocalRunner, SlurmRunner
    >>>
    >>> # For local execution
    >>> runner = LocalRunner(vasp_command='vasp_std')
    >>>
    >>> # For SLURM cluster
    >>> runner = SlurmRunner(partition='compute', nodes=2)
"""

from .base import Runner, JobState, JobStatus
from .local import LocalRunner
from .mock import MockRunner, MockResults
from .slurm import SlurmRunner
from .kubernetes import KubernetesRunner
from .interactive import InteractiveRunner, InteractiveResults
from .socket_io import SocketServer, SocketClient, SocketConfig

__all__ = [
    'Runner',
    'JobState',
    'JobStatus',
    'LocalRunner',
    'MockRunner',
    'MockResults',
    'SlurmRunner',
    'KubernetesRunner',
    'InteractiveRunner',
    'InteractiveResults',
    'SocketServer',
    'SocketClient',
    'SocketConfig',
]

# Optional runners imported on demand
def get_slurm_runner():
    """Get SlurmRunner class (avoids import if not needed)."""
    from .slurm import SlurmRunner
    return SlurmRunner

def get_kubernetes_runner():
    """Get KubernetesRunner class (requires kubernetes package)."""
    from .kubernetes import KubernetesRunner
    return KubernetesRunner

def get_interactive_runner():
    """Get InteractiveRunner class for persistent VASP sessions."""
    from .interactive import InteractiveRunner
    return InteractiveRunner
