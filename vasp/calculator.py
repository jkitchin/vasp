"""Main VASP calculator class.

This module provides the Vasp class, which is the primary interface
for running VASP calculations through ASE. It combines:
- ASE FileIOCalculator base class
- Mixin classes for specific functionality
- Pluggable runner system for different execution backends
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from .exceptions import (
    VaspException,
    VaspSubmitted,
    VaspQueued,
    VaspRunning,
    VaspNotFinished,
    VaspNotConverged,
    VaspSetupError,
)
from .mixins import IOMixin, ElectronicMixin, AnalysisMixin, DynamicsMixin
from .runners import Runner, LocalRunner, JobState, JobStatus

if TYPE_CHECKING:
    from ase import Atoms

# Configure logging
log = logging.getLogger('vasp')


# Exchange-correlation functional settings
XC_DEFAULTS: dict[str, dict[str, Any]] = {
    'lda': {'pp': 'LDA'},
    'pbe': {'gga': 'PE'},
    'revpbe': {'gga': 'RE'},
    'rpbe': {'gga': 'RP'},
    'pbesol': {'gga': 'PS'},
    'am05': {'gga': 'AM'},
    'pbe0': {'gga': 'PE', 'lhfcalc': True, 'aexx': 0.25},
    'hse03': {'gga': 'PE', 'lhfcalc': True, 'hfscreen': 0.3},
    'hse06': {'gga': 'PE', 'lhfcalc': True, 'hfscreen': 0.2},
    'b3lyp': {'gga': 'B3', 'lhfcalc': True, 'aexx': 0.2, 'aggax': 0.72, 'aggac': 0.81, 'aldac': 0.19},
    'scan': {'metagga': 'SCAN'},
    'tpss': {'metagga': 'TPSS'},
    'optpbe-vdw': {'gga': 'OR', 'luse_vdw': True, 'aggac': 0.0},
    'optb88-vdw': {'gga': 'BO', 'luse_vdw': True, 'aggac': 0.0, 'param1': 1.1/6.0, 'param2': 0.22},
    'optb86b-vdw': {'gga': 'MK', 'luse_vdw': True, 'aggac': 0.0, 'param1': 0.1234, 'param2': 1.0},
    'vdw-df2': {'gga': 'ML', 'luse_vdw': True, 'aggac': 0.0, 'zab_vdw': -1.8867},
    'beef-vdw': {'gga': 'BF', 'luse_vdw': True, 'zab_vdw': -1.8867, 'lbeefens': True},
}


@dataclass
class CalculationResult:
    """Container for workflow-friendly calculation results.

    Use with `calc.run_async()` for exception-free status checking.

    Attributes:
        state: Current job state.
        energy: Total energy in eV (if complete).
        forces: Forces array (if complete).
        stress: Stress tensor (if complete).
        jobid: Job identifier (if submitted).
        error: Error message (if failed).
    """
    state: JobState
    energy: float | None = None
    forces: np.ndarray | None = None
    stress: np.ndarray | None = None
    jobid: str | None = None
    error: str | None = None


class Vasp(Calculator, IOMixin, ElectronicMixin, AnalysisMixin, DynamicsMixin):
    """ASE calculator interface for VASP.

    This calculator provides a complete interface to VASP through ASE,
    with support for:
    - All standard VASP parameters as keyword arguments
    - Automatic input file generation
    - Result parsing from output files
    - Pluggable execution backends (local, SLURM, Kubernetes)
    - Non-blocking async execution with exception-based signaling

    Args:
        label: Calculation directory path (default: 'vasp').
        atoms: ASE Atoms object for the calculation.
        runner: Execution backend (default: LocalRunner).
        xc: Exchange-correlation functional (e.g., 'PBE', 'HSE06').
        pp: Pseudopotential type (e.g., 'PBE', 'LDA').
        kpts: K-point grid as (nx, ny, nz) tuple.
        gamma: Use Gamma-centered k-point grid.
        setups: Dict of special POTCAR setups per element.
        **kwargs: Any VASP INCAR parameters.

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
        >>>
        >>> try:
        ...     energy = calc.potential_energy
        ... except VaspSubmitted as e:
        ...     print(f"Job submitted: {e.jobid}")
        ... except VaspQueued:
        ...     print("Waiting in queue...")
    """

    name = 'vasp'
    implemented_properties = ['energy', 'forces', 'stress', 'magmom', 'magmoms']
    default_parameters: dict[str, Any] = {
        'xc': 'PBE',
        'pp': 'PBE',
        'kpts': (1, 1, 1),
        'gamma': False,
        'ismear': 1,
        'sigma': 0.1,
        'lwave': False,
        'lcharg': False,
    }

    def __init__(
        self,
        label: str = 'vasp',
        atoms: Atoms | None = None,
        runner: Runner | None = None,
        **kwargs
    ):
        # Initialize parent Calculator first
        Calculator.__init__(self, atoms=atoms)

        # Set up runner
        if runner is None:
            runner = LocalRunner()
        self.runner = runner

        # Process directory from label (after parent init to avoid being overwritten)
        self.directory = os.path.abspath(label)

        # Initialize parameters
        self.parameters: dict[str, Any] = {}

        # Merge default parameters
        for key, val in self.default_parameters.items():
            self.parameters[key] = val

        # Process kwargs
        self._process_parameters(kwargs)

        # Store atoms and set up sorting
        if atoms is not None:
            self.atoms = atoms
            self._setup_sorting(atoms)
        else:
            self.sort = []
            self.resort = []

        # Results storage
        self.results: dict[str, Any] = {}

    def _process_parameters(self, kwargs: dict) -> None:
        """Process and validate input parameters."""
        for key, val in kwargs.items():
            key_lower = key.lower()

            # Handle XC functional
            if key_lower == 'xc':
                self.parameters['xc'] = val
                xc_lower = val.lower()
                if xc_lower in XC_DEFAULTS:
                    for xc_key, xc_val in XC_DEFAULTS[xc_lower].items():
                        if xc_key not in kwargs:
                            self.parameters[xc_key] = xc_val
            else:
                self.parameters[key_lower] = val

    def _setup_sorting(self, atoms: Atoms) -> None:
        """Set up atom sorting by chemical symbol.

        VASP requires atoms to be grouped by species. This creates
        sort/resort indices to handle the reordering.
        """
        symbols = atoms.get_chemical_symbols()

        # Group by symbol
        symbol_order = []
        for s in symbols:
            if s not in symbol_order:
                symbol_order.append(s)

        # Create sorting indices
        self.sort = []
        for symbol in symbol_order:
            for i, s in enumerate(symbols):
                if s == symbol:
                    self.sort.append(i)

        # Create reverse sorting
        self.resort = [0] * len(self.sort)
        for i, j in enumerate(self.sort):
            self.resort[j] = i

    def set(self, **kwargs) -> dict:
        """Set calculator parameters.

        Args:
            **kwargs: Parameters to update.

        Returns:
            Dict of changed parameters.
        """
        changed = {}

        for key, val in kwargs.items():
            key_lower = key.lower()
            old_val = self.parameters.get(key_lower)

            if old_val != val:
                changed[key_lower] = val
                self.parameters[key_lower] = val

                # Handle XC changes
                if key_lower == 'xc':
                    xc_lower = val.lower()
                    if xc_lower in XC_DEFAULTS:
                        for xc_key, xc_val in XC_DEFAULTS[xc_lower].items():
                            self.parameters[xc_key] = xc_val

        if changed:
            self.results = {}

        return changed

    def calculate(
        self,
        atoms: Atoms | None = None,
        properties: list[str] | None = None,
        system_changes: list[str] | None = None,
    ) -> None:
        """Run VASP calculation.

        This method triggers the VASP calculation. Depending on the
        runner configuration, it may:
        - Run VASP and wait for completion (LocalRunner, blocking)
        - Submit job and raise exception (SLURM, K8s runners)

        Args:
            atoms: ASE Atoms object (uses self.atoms if None).
            properties: Properties to calculate.
            system_changes: What changed since last calculation.

        Raises:
            VaspSubmitted: Job was just submitted.
            VaspQueued: Job is in queue.
            VaspRunning: Job is currently running.
            VaspNotConverged: Calculation failed to converge.
        """
        if properties is None:
            properties = self.implemented_properties

        if atoms is not None:
            self.atoms = atoms
            self._setup_sorting(atoms)

        # Check current status
        status = self.runner.status(self.directory)
        log.debug(f"Current status: {status.state}")

        if status.state == JobState.COMPLETE:
            self.read_results()
            return

        if status.state == JobState.QUEUED:
            raise VaspQueued(message=status.message, jobid=status.jobid)

        if status.state == JobState.RUNNING:
            raise VaspRunning(message=status.message, jobid=status.jobid)

        if status.state == JobState.FAILED:
            raise VaspNotConverged(status.message or "Calculation failed")

        # Not started - need to run
        Calculator.calculate(self, atoms, properties, system_changes)

        # Write input files
        self.write_input(atoms, properties, system_changes)

        # Run calculation
        result = self.runner.run(self.directory)

        if result.state == JobState.COMPLETE:
            self.read_results()

    def update(self) -> None:
        """Ensure calculation results are current.

        Checks status and reads results if complete.
        Raises exceptions for non-complete states.
        """
        status = self.runner.status(self.directory)

        if status.state == JobState.COMPLETE:
            if not self.results:
                self.read_results()
            return

        if status.state == JobState.QUEUED:
            raise VaspQueued(message=status.message, jobid=status.jobid)

        if status.state == JobState.RUNNING:
            raise VaspRunning(message=status.message, jobid=status.jobid)

        if status.state == JobState.FAILED:
            raise VaspNotConverged(status.message or "Calculation failed")

        if status.state == JobState.NOT_STARTED:
            raise VaspNotFinished("Calculation not started")

    # =========================================================================
    # ASE Calculator Interface Properties
    # =========================================================================

    @property
    def potential_energy(self) -> float:
        """Get potential energy in eV."""
        self.update()
        return self.results['energy']

    @property
    def forces(self) -> np.ndarray:
        """Get forces in eV/Angstrom."""
        self.update()
        return self.results['forces']

    @property
    def stress(self) -> np.ndarray:
        """Get stress tensor in eV/Angstrom^3."""
        self.update()
        return self.results['stress']

    def get_potential_energy(self, atoms: Atoms | None = None, force_consistent: bool = False) -> float:
        """Get potential energy.

        Args:
            atoms: Atoms object (triggers calculation if different).
            force_consistent: If True, return energy consistent with forces.

        Returns:
            Total energy in eV.
        """
        if atoms is not None:
            self.calculate(atoms=atoms, properties=['energy'])
        else:
            self.update()

        if force_consistent and 'free_energy' in self.results:
            return self.results['free_energy']
        return self.results['energy']

    def get_forces(self, atoms: Atoms | None = None) -> np.ndarray:
        """Get forces on atoms.

        Args:
            atoms: Atoms object (triggers calculation if different).

        Returns:
            Forces array of shape (natoms, 3) in eV/Angstrom.
        """
        if atoms is not None:
            self.calculate(atoms=atoms, properties=['forces'])
        else:
            self.update()

        return self.results['forces']

    def get_stress(self, atoms: Atoms | None = None) -> np.ndarray:
        """Get stress tensor.

        Args:
            atoms: Atoms object (triggers calculation if different).

        Returns:
            Stress in Voigt notation (6,) in eV/Angstrom^3.
        """
        if atoms is not None:
            self.calculate(atoms=atoms, properties=['stress'])
        else:
            self.update()

        return self.results['stress']

    # =========================================================================
    # Workflow-Friendly Methods
    # =========================================================================

    def submit(self) -> str | None:
        """Submit calculation without blocking.

        Writes input files and submits to runner.
        Returns job ID if applicable.

        Returns:
            Job ID string, or None for local runner.
        """
        self.write_input(self.atoms)

        try:
            self.runner.run(self.directory)
            return None
        except VaspSubmitted as e:
            return e.jobid
        except (VaspQueued, VaspRunning) as e:
            return e.jobid

    def poll(self) -> JobStatus:
        """Check calculation status without triggering anything.

        Returns:
            JobStatus with current state.
        """
        return self.runner.status(self.directory)

    def is_complete(self) -> bool:
        """Check if calculation is complete.

        Returns:
            True if calculation finished successfully.
        """
        return self.poll().state == JobState.COMPLETE

    def run_async(self) -> CalculationResult:
        """Run calculation with result object instead of exceptions.

        Useful for workflow tools that prefer return values.

        Returns:
            CalculationResult with state and results.

        Example:
            >>> result = calc.run_async()
            >>> if result.state == JobState.COMPLETE:
            ...     print(result.energy)
        """
        status = self.runner.status(self.directory)

        if status.state == JobState.COMPLETE:
            self.read_results()
            return CalculationResult(
                state=JobState.COMPLETE,
                energy=self.results.get('energy'),
                forces=self.results.get('forces'),
                stress=self.results.get('stress'),
            )

        if status.state == JobState.NOT_STARTED:
            self.write_input(self.atoms)
            try:
                new_status = self.runner.run(self.directory)
                # If completed immediately, read results
                if new_status.state == JobState.COMPLETE:
                    self.read_results()
                    return CalculationResult(
                        state=JobState.COMPLETE,
                        energy=self.results.get('energy'),
                        forces=self.results.get('forces'),
                        stress=self.results.get('stress'),
                    )
                return CalculationResult(
                    state=new_status.state,
                    jobid=new_status.jobid,
                )
            except VaspSubmitted as e:
                return CalculationResult(
                    state=JobState.SUBMITTED,
                    jobid=e.jobid,
                )
            except (VaspQueued, VaspRunning) as e:
                return CalculationResult(
                    state=JobState.RUNNING,
                    jobid=e.jobid,
                )

        return CalculationResult(
            state=status.state,
            jobid=status.jobid,
            error=status.message,
        )

    def cancel(self) -> bool:
        """Cancel running calculation.

        Returns:
            True if cancellation successful.
        """
        return self.runner.cancel(self.directory)

    # =========================================================================
    # Serialization
    # =========================================================================

    def __getstate__(self) -> dict:
        """Prepare for pickling (for workflow tools)."""
        state = self.__dict__.copy()
        # Don't pickle runner - recreate on unpickle
        state['_runner_type'] = type(self.runner).__name__
        state['runner'] = None
        return state

    def __setstate__(self, state: dict) -> None:
        """Restore from pickle."""
        self.__dict__.update(state)
        # Recreate default runner
        self.runner = LocalRunner()

    def __repr__(self) -> str:
        return f"Vasp('{self.directory}', xc='{self.parameters.get('xc', 'PBE')}')"
