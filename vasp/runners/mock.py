"""Mock runner for testing VASP calculations.

This runner simulates VASP execution by creating fake output files
from predefined templates or fixture data. It's designed for:
- Unit testing calculator logic without running VASP
- Integration testing workflow pipelines
- Development when VASP is not available

Example:
    >>> from vasp.runners import MockRunner
    >>>
    >>> # Simple mock that always succeeds with fixed energy
    >>> runner = MockRunner(energy=-10.5, forces=[[0, 0, 0], [0, 0, 0.1]])
    >>>
    >>> calc = Vasp('test_calc', runner=runner, atoms=atoms)
    >>> energy = calc.potential_energy  # Returns -10.5
    >>>
    >>> # Mock with custom state sequence for testing async workflows
    >>> runner = MockRunner(
    ...     state_sequence=[JobState.QUEUED, JobState.RUNNING, JobState.COMPLETE],
    ...     energy=-10.5,
    ... )
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

from ..exceptions import VaspQueued, VaspRunning, VaspSubmitted
from .base import JobState, JobStatus, Runner

if TYPE_CHECKING:
    pass


@dataclass
class MockResults:
    """Container for mock calculation results.

    Attributes:
        energy: Total energy in eV.
        forces: Forces array of shape (natoms, 3) in eV/Angstrom.
        stress: Stress tensor as 6-element Voigt array in kBar.
        magmom: Total magnetic moment in Bohr magnetons.
        magmoms: Per-atom magnetic moments.
        fermi_level: Fermi energy in eV.
        eigenvalues: Eigenvalue array (nspin, nkpts, nbands).
        occupations: Occupation numbers (nspin, nkpts, nbands).
        ibz_kpts: IBZ k-points array (nkpts, 3).
        kpt_weights: K-point weights (nkpts,).
        converged: Whether calculation converged.
        metadata: Additional data.
    """
    energy: float = -10.0
    forces: list[list[float]] | np.ndarray | None = None
    stress: list[float] | np.ndarray | None = None
    magmom: float = 0.0
    magmoms: list[float] | np.ndarray | None = None
    fermi_level: float = -5.0
    eigenvalues: np.ndarray | None = None
    occupations: np.ndarray | None = None
    ibz_kpts: np.ndarray | None = None
    kpt_weights: np.ndarray | None = None
    converged: bool = True
    # Additional properties
    metadata: dict[str, Any] = field(default_factory=dict)


# Global counter for job IDs
_mock_job_counter = 0


class MockRunner(Runner):
    """Mock runner for testing without VASP.

    Simulates VASP execution by writing fake output files.
    Useful for testing calculator logic, workflow integration,
    and development when VASP is not available.

    Args:
        results: MockResults object with calculation outputs.
        energy: Shortcut to set results.energy.
        forces: Shortcut to set results.forces.
        state_sequence: List of JobStates to cycle through on each
            status() call. Useful for testing async workflows.
        delay_calls: Number of status() calls before returning COMPLETE.
            Alternative to state_sequence for simple cases.
        fail_on_run: If True, immediately fail when run() is called.
        error_message: Error message if fail_on_run is True.
        write_outputs: If True, write mock OUTCAR/vasprun.xml files.

    Example:
        >>> # Simple successful calculation
        >>> runner = MockRunner(energy=-15.3)
        >>>
        >>> # Simulate job that takes 3 status checks to complete
        >>> runner = MockRunner(
        ...     energy=-15.3,
        ...     state_sequence=[
        ...         JobState.QUEUED,
        ...         JobState.RUNNING,
        ...         JobState.RUNNING,
        ...         JobState.COMPLETE,
        ...     ]
        ... )
        >>>
        >>> # Simulate failed calculation
        >>> runner = MockRunner(fail_on_run=True, error_message="ZBRENT error")
    """

    def __init__(
        self,
        results: MockResults | None = None,
        energy: float | None = None,
        forces: list | np.ndarray | None = None,
        stress: list | np.ndarray | None = None,
        state_sequence: list[JobState] | None = None,
        delay_calls: int = 0,
        delay: float = 0,
        fail_on_run: bool = False,
        error_message: str = "Mock error",
        write_outputs: bool = True,
    ):
        # Build results from shortcuts
        if results is None:
            results = MockResults()
        if energy is not None:
            results.energy = energy
        if forces is not None:
            results.forces = forces
        if stress is not None:
            results.stress = stress

        self.results = results
        self.state_sequence = state_sequence
        self.delay_calls = delay_calls
        self.fail_on_run = fail_on_run
        self.error_message = error_message
        self.write_outputs = write_outputs
        self.delay = delay

        # Track state for each directory
        self._call_counts: dict[str, int] = {}
        self._states: dict[str, JobState] = {}
        self._job_ids: dict[str, str] = {}

        # Get unique job ID for this runner
        global _mock_job_counter
        _mock_job_counter += 1
        self._runner_id = _mock_job_counter

    def _get_job_id(self, directory: str) -> str:
        """Get consistent job ID for a directory."""
        if directory not in self._job_ids:
            self._job_ids[directory] = f"mock-job-{self._runner_id}"
        return self._job_ids[directory]

    def run(self, directory: str) -> JobStatus:
        """Simulate running VASP.

        Creates mock output files and tracks job state.
        """
        import time

        if self.fail_on_run:
            self._states[directory] = JobState.FAILED
            return JobStatus(JobState.FAILED, message=self.error_message)

        # Apply delay if specified
        if self.delay > 0:
            time.sleep(self.delay)

        # Initialize tracking - start at 1 so status() returns next state in sequence
        self._call_counts[directory] = 1
        job_id = self._get_job_id(directory)

        if self.state_sequence:
            # Use first state in sequence
            state = self.state_sequence[0]
            self._states[directory] = state

            if state == JobState.SUBMITTED:
                raise VaspSubmitted(jobid=job_id)
            elif state == JobState.QUEUED:
                raise VaspQueued(jobid=job_id)
            elif state == JobState.RUNNING:
                raise VaspRunning(jobid=job_id)
        elif self.delay_calls > 0:
            # Start as running
            self._states[directory] = JobState.RUNNING
            raise VaspRunning(jobid=job_id)
        else:
            # Immediate completion
            self._states[directory] = JobState.COMPLETE
            if self.write_outputs:
                self._write_mock_outputs(directory)
            return JobStatus(JobState.COMPLETE)

        return JobStatus(self._states[directory])

    def status(self, directory: str) -> JobStatus:
        """Check mock job status.

        Cycles through state_sequence or delay_calls to simulate
        job progression.
        """
        if directory not in self._states:
            return JobStatus(JobState.NOT_STARTED)

        count = self._call_counts.get(directory, 0)
        self._call_counts[directory] = count + 1
        job_id = self._get_job_id(directory)

        if self.state_sequence:
            # Cycle through state sequence
            idx = min(count, len(self.state_sequence) - 1)
            state = self.state_sequence[idx]
            self._states[directory] = state

            if state == JobState.COMPLETE and self.write_outputs:
                self._write_mock_outputs(directory)

            return JobStatus(state, jobid=job_id)

        elif self.delay_calls > 0:
            if count >= self.delay_calls:
                self._states[directory] = JobState.COMPLETE
                if self.write_outputs:
                    self._write_mock_outputs(directory)
                return JobStatus(JobState.COMPLETE)
            else:
                return JobStatus(JobState.RUNNING, jobid=job_id)

        return JobStatus(self._states[directory])

    def cancel(self, directory: str) -> bool:
        """Cancel mock job."""
        if directory in self._states:
            self._states[directory] = JobState.FAILED
        return True

    def reset(self) -> None:
        """Reset all state tracking (useful between tests)."""
        self._call_counts.clear()
        self._states.clear()
        self._job_ids.clear()

    def _write_mock_outputs(self, directory: str) -> None:
        """Write mock VASP output files."""
        os.makedirs(directory, exist_ok=True)

        # Write mock OUTCAR
        self._write_mock_outcar(directory)

        # Write mock vasprun.xml
        self._write_mock_vasprun(directory)

        # Write mock CONTCAR (copy of POSCAR if exists)
        poscar = os.path.join(directory, 'POSCAR')
        contcar = os.path.join(directory, 'CONTCAR')
        if os.path.exists(poscar) and not os.path.exists(contcar):
            with open(poscar) as f:
                content = f.read()
            with open(contcar, 'w') as f:
                f.write(content)

    def _write_mock_outcar(self, directory: str) -> None:
        """Write mock OUTCAR file."""
        r = self.results
        forces_str = ""
        if r.forces is not None:
            forces = np.asarray(r.forces)
            forces_str = " POSITION                                       TOTAL-FORCE (eV/Angst)\n"
            forces_str += " -----------------------------------------------------------------------------------\n"
            for i, f in enumerate(forces):
                forces_str += f"      0.00000      0.00000      {i:.5f}         {f[0]:.6f}      {f[1]:.6f}      {f[2]:.6f}\n"
            forces_str += " -----------------------------------------------------------------------------------\n"

        stress_str = ""
        if r.stress is not None:
            s = np.asarray(r.stress)
            stress_str = f"""  in kB  {s[0]:.2f}  {s[1]:.2f}  {s[2]:.2f}  {s[3]:.2f}  {s[4]:.2f}  {s[5]:.2f}\n"""

        # Magnetic moment line
        magmom_str = f" number of electron       8.0000000 magnetization       {r.magmom:.7f}\n"

        # Per-atom magnetic moments
        magmoms_str = ""
        if r.magmoms is not None:
            magmoms = np.asarray(r.magmoms)
            magmoms_str = " magnetization (x)\n\n"
            for i, m in enumerate(magmoms):
                magmoms_str += f"     {i + 1}     {m:.5f}\n"

        content = f"""
 vasp.6.3.0 (mock output for testing)
 executed on             LinuxIFC date 2024.01.01  00:00:00

 POTCAR:    PAW_PBE Si 05Jan2001

 POSCAR: Mock structure

 Startparameter for this run:
   ENCUT  =  400.0 eV

 E-fermi :   {r.fermi_level:.4f}     XC(G=0):  -0.0000     alpha+bet : -0.0000

{forces_str}
{stress_str}
  free  energy   TOTEN  =      {r.energy:.8f} eV

  energy  without entropy=      {r.energy:.8f}  energy(sigma->0) =      {r.energy:.8f}

{magmom_str}
{magmoms_str}

 # of diffusion steps           =      0

 General timing and accounting informance for this job:

       LOOP:  cpu time  123.45: real time  123.45
"""
        outcar_path = os.path.join(directory, 'OUTCAR')
        with open(outcar_path, 'w') as f:
            f.write(content)

    def _write_mock_vasprun(self, directory: str) -> None:
        """Write mock vasprun.xml file."""
        r = self.results
        forces_xml = ""
        if r.forces is not None:
            forces = np.asarray(r.forces)
            forces_xml = '<varray name="forces">\n'
            for f in forces:
                forces_xml += f'  <v>  {f[0]:.8f}  {f[1]:.8f}  {f[2]:.8f} </v>\n'
            forces_xml += '</varray>\n'

        stress_xml = ""
        if r.stress is not None:
            s = np.asarray(r.stress)
            # Convert Voigt to 3x3
            stress_xml = '<varray name="stress">\n'
            stress_xml += f'  <v>  {s[0]:.8f}  {s[5]:.8f}  {s[4]:.8f} </v>\n'
            stress_xml += f'  <v>  {s[5]:.8f}  {s[1]:.8f}  {s[3]:.8f} </v>\n'
            stress_xml += f'  <v>  {s[4]:.8f}  {s[3]:.8f}  {s[2]:.8f} </v>\n'
            stress_xml += '</varray>\n'

        # Build eigenvalues XML if provided
        eigenvalues_xml = ""
        if r.eigenvalues is not None:
            eigenvalues = np.asarray(r.eigenvalues)
            occupations = np.asarray(r.occupations) if r.occupations is not None else np.ones_like(eigenvalues)
            eigenvalues_xml = '<eigenvalues>\n<array>\n<set>\n'
            for spin_idx in range(eigenvalues.shape[0]):
                eigenvalues_xml += f'  <set comment="spin {spin_idx + 1}">\n'
                for kpt_idx in range(eigenvalues.shape[1]):
                    eigenvalues_xml += f'    <set comment="kpoint {kpt_idx + 1}">\n'
                    for band_idx in range(eigenvalues.shape[2]):
                        eig = eigenvalues[spin_idx, kpt_idx, band_idx]
                        occ = occupations[spin_idx, kpt_idx, band_idx]
                        eigenvalues_xml += f'      <r>  {eig:.8f}  {occ:.8f}</r>\n'
                    eigenvalues_xml += '    </set>\n'
                eigenvalues_xml += '  </set>\n'
            eigenvalues_xml += '</set>\n</array>\n</eigenvalues>\n'

        # Build k-points XML if provided
        kpoints_xml = ""
        if r.ibz_kpts is not None:
            ibz_kpts = np.asarray(r.ibz_kpts)
            kpt_weights = np.asarray(r.kpt_weights) if r.kpt_weights is not None else np.ones(len(ibz_kpts)) / len(ibz_kpts)
            kpoints_xml = '<kpoints>\n<varray name="kpointlist">\n'
            for kpt in ibz_kpts:
                kpoints_xml += f'  <v>  {kpt[0]:.8f}  {kpt[1]:.8f}  {kpt[2]:.8f}</v>\n'
            kpoints_xml += '</varray>\n<varray name="weights">\n'
            for w in kpt_weights:
                kpoints_xml += f'  <v>  {w:.8f}</v>\n'
            kpoints_xml += '</varray>\n</kpoints>\n'

        content = f"""<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <generator>
  <i name="program" type="string">vasp </i>
  <i name="version" type="string">6.3.0</i>
 </generator>
 {kpoints_xml}
 <calculation>
  <energy>
   <i name="e_fr_energy">    {r.energy:.8f} </i>
   <i name="e_wo_entrp">    {r.energy:.8f} </i>
   <i name="e_0_energy">    {r.energy:.8f} </i>
  </energy>
  {forces_xml}
  {stress_xml}
  {eigenvalues_xml}
  <dos>
   <i name="efermi">    {r.fermi_level:.8f} </i>
  </dos>
 </calculation>
</modeling>
"""
        vasprun_path = os.path.join(directory, 'vasprun.xml')
        with open(vasprun_path, 'w') as f:
            f.write(content)

    def __repr__(self) -> str:
        return (
            f"MockRunner(energy={self.results.energy}, "
            f"fail_on_run={self.fail_on_run})"
        )
