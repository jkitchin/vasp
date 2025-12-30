"""VASP calculator interface for ASE.

This package provides a modern Python interface to VASP (Vienna Ab initio
Simulation Package) through ASE (Atomic Simulation Environment).

Features:
- Full VASP parameter support via keyword arguments
- Automatic input file generation
- Result parsing from output files
- Pluggable execution backends (local, SLURM, Kubernetes)
- Non-blocking async execution with exception-based signaling
- Workflow tool integration (Prefect, Dask, Parsl, Covalent)
- Parameter presets for common calculation types
- quacc-style recipes for automated workflows

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

Using parameter presets:
    >>> from vasp.parameters import get_vdw_params, get_ldau_params
    >>>
    >>> # Add D3-BJ dispersion correction
    >>> vdw = get_vdw_params('d3bj')
    >>> calc = Vasp(..., **vdw)
    >>>
    >>> # DFT+U for Fe oxide
    >>> ldau = get_ldau_params(['Fe', 'O'], {'Fe': 4.0})
    >>> calc = Vasp(..., **ldau)

Using recipes:
    >>> from vasp.recipes import relax_job, static_job, phonon_flow
    >>>
    >>> result = relax_job(atoms, relax_cell=True)
    >>> phonons = phonon_flow(result.atoms, supercell_matrix=(2, 2, 2))
"""

__version__ = '2.0.0'
__author__ = 'John Kitchin'

# Main calculator class
from .calculator import CalculationResult, Vasp

# Exceptions
from .exceptions import (
    VaspEmptyOutput,
    VaspError,
    VaspException,
    VaspNotConverged,
    VaspNotFinished,
    VaspQueued,
    VaspRunning,
    VaspSetupError,
    VaspSubmitted,
    VaspWarning,
)

# Runners
from .runners import (
    JobState,
    JobStatus,
    LocalRunner,
    MockRunner,
    Runner,
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
