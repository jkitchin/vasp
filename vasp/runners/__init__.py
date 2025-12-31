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

from .base import JobState, JobStatus, Runner
from .interactive import InteractiveResults, InteractiveRunner
from .kubernetes import KubernetesRunner
from .local import LocalRunner, get_optimal_nprocs
from .mock import MockResults, MockRunner
from .slurm import SlurmRunner
from .socket_io import SocketClient, SocketConfig, SocketServer

__all__ = [
    "Runner",
    "JobState",
    "JobStatus",
    "LocalRunner",
    "get_optimal_nprocs",
    "MockRunner",
    "MockResults",
    "SlurmRunner",
    "KubernetesRunner",
    "InteractiveRunner",
    "InteractiveResults",
    "SocketServer",
    "SocketClient",
    "SocketConfig",
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
