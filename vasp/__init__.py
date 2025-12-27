"""VASP calculator interface for ASE.

This package provides a modern Python interface to VASP (Vienna Ab initio
Simulation Package) through ASE (Atomic Simulation Environment).

Features:
- Full VASP parameter support via keyword arguments
- Automatic input file generation
- Result parsing from output files
- Pluggable execution backends (local, SLURM, Kubernetes)
- Non-blocking async execution with exception-based signaling
- Workflow tool integration (Prefect, Dask, etc.)

Example:
    >>> from ase.build import bulk
    >>> from vasp import Vasp
    >>>
    >>> atoms = bulk('Si')
    >>> calc = Vasp(
    ...     'si_calc',
    ...     atoms=atoms,
    ...     xc='PBE',
    ...     encut=400,
    ...     kpts=(8, 8, 8),
    ... )
    >>>
    >>> energy = calc.potential_energy
    >>> print(f"Energy: {energy:.3f} eV")

Non-blocking usage:
    >>> from vasp.exceptions import VaspSubmitted, VaspQueued
    >>> from vasp.runners import SlurmRunner
    >>>
    >>> runner = SlurmRunner(partition='compute', nodes=2)
    >>> calc = Vasp('my_calc', atoms=atoms, runner=runner)
    >>>
    >>> try:
    ...     energy = calc.potential_energy
    ... except VaspSubmitted as e:
    ...     print(f"Job submitted: {e.jobid}")
"""

__version__ = '1.0.0'
__author__ = 'John Kitchin'

# Main calculator class
from .calculator import Vasp, CalculationResult

# Exceptions
from .exceptions import (
    VaspException,
    VaspSubmitted,
    VaspQueued,
    VaspRunning,
    VaspNotFinished,
    VaspNotConverged,
    VaspError,
    VaspEmptyOutput,
    VaspSetupError,
    VaspWarning,
)

# Runners
from .runners import (
    Runner,
    JobState,
    JobStatus,
    LocalRunner,
    MockRunner,
)

__all__ = [
    # Version
    '__version__',
    # Main class
    'Vasp',
    'CalculationResult',
    # Exceptions
    'VaspException',
    'VaspSubmitted',
    'VaspQueued',
    'VaspRunning',
    'VaspNotFinished',
    'VaspNotConverged',
    'VaspError',
    'VaspEmptyOutput',
    'VaspSetupError',
    'VaspWarning',
    # Runners
    'Runner',
    'JobState',
    'JobStatus',
    'LocalRunner',
    'MockRunner',
]


def get_slurm_runner():
    """Get SlurmRunner class.

    Returns:
        SlurmRunner class.

    Example:
        >>> SlurmRunner = get_slurm_runner()
        >>> runner = SlurmRunner(partition='compute')
    """
    from .runners.slurm import SlurmRunner
    return SlurmRunner


def get_kubernetes_runner():
    """Get KubernetesRunner class.

    Requires: pip install kubernetes

    Returns:
        KubernetesRunner class.

    Example:
        >>> KubernetesRunner = get_kubernetes_runner()
        >>> runner = KubernetesRunner(namespace='vasp-jobs')
    """
    from .runners.kubernetes import KubernetesRunner
    return KubernetesRunner
