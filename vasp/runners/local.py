"""Local runner for VASP execution.

Runs VASP directly on the local machine, either blocking (default)
or in background mode for non-blocking operation.
"""

from __future__ import annotations

import os
import signal
import subprocess
from typing import TYPE_CHECKING

from ..exceptions import VaspRunning, VaspSetupError
from .base import JobState, JobStatus, Runner

if TYPE_CHECKING:
    pass


class LocalRunner(Runner):
    """Run VASP on the local machine.

    By default runs synchronously (blocking until complete).
    Use background=True for non-blocking execution.

    Args:
        vasp_command: Command to run VASP. Defaults to $VASP_COMMAND
            environment variable, or 'vasp_std' if not set.
        mpi_command: MPI launcher command, e.g., 'mpirun -np 4'.
            If None, runs serial VASP.
        background: If True, run in background and return immediately.
        working_dir_strategy: How to handle working directory.
            'chdir' = change to directory before running.
            'arg' = pass directory as argument (not common for VASP).

    Example:
        >>> runner = LocalRunner(vasp_command='vasp_std')
        >>> calc = Vasp('my_calc', runner=runner, ...)
        >>> energy = calc.potential_energy  # Blocks until done

        >>> # Non-blocking
        >>> runner = LocalRunner(vasp_command='vasp_std', background=True)
        >>> calc = Vasp('my_calc', runner=runner, ...)
        >>> try:
        ...     energy = calc.potential_energy
        ... except VaspRunning:
        ...     print("Still running...")
    """

    def __init__(
        self,
        vasp_command: str | None = None,
        mpi_command: str | None = None,
        background: bool = False,
    ):
        self.vasp_command = vasp_command or os.environ.get('VASP_COMMAND', 'vasp_std')
        self.mpi_command = mpi_command
        self.background = background

    def run(self, directory: str) -> JobStatus:
        """Run VASP in the specified directory.

        Args:
            directory: Path to calculation directory with input files.

        Returns:
            JobStatus indicating outcome.

        Raises:
            VaspRunning: If background=True and job is still running.
            VaspSetupError: If required input files are missing.
        """
        # Check if already complete
        current = self.status(directory)
        if current.state == JobState.COMPLETE:
            return current

        # Check if already running
        if current.state == JobState.RUNNING:
            raise VaspRunning(jobid=current.jobid)

        # Verify input files exist
        self._verify_inputs(directory)

        # Build command
        cmd = self._build_command()

        if self.background:
            pid = self._run_background(directory, cmd)
            raise VaspRunning(
                message=f"Started background process {pid}",
                jobid=str(pid)
            )
        else:
            return self._run_blocking(directory, cmd)

    def status(self, directory: str) -> JobStatus:
        """Check status of calculation."""
        # Check if process is running
        pid = self._read_pid(directory)
        if pid and self._is_process_running(pid):
            return JobStatus(JobState.RUNNING, jobid=str(pid))

        # Check OUTCAR for completion
        if self._check_outcar_complete(directory):
            return JobStatus(JobState.COMPLETE)

        # Check for errors
        error = self._check_outcar_error(directory)
        if error:
            return JobStatus(JobState.FAILED, message=error)

        # Check if OUTCAR exists but incomplete
        outcar = os.path.join(directory, 'OUTCAR')
        if os.path.exists(outcar):
            # OUTCAR exists but not complete and no process
            return JobStatus(
                JobState.FAILED,
                message="OUTCAR incomplete and no running process"
            )

        return JobStatus(JobState.NOT_STARTED)

    def cancel(self, directory: str) -> bool:
        """Kill running VASP process."""
        pid = self._read_pid(directory)
        if pid and self._is_process_running(pid):
            try:
                os.kill(pid, signal.SIGTERM)
                return True
            except OSError:
                return False
        return True

    def _build_command(self) -> str:
        """Build the full command string."""
        if self.mpi_command:
            return f"{self.mpi_command} {self.vasp_command}"
        return self.vasp_command

    def _run_blocking(self, directory: str, cmd: str) -> JobStatus:
        """Run VASP and wait for completion."""
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=directory,
                capture_output=True,
                text=True,
            )

            if result.returncode == 0 and self._check_outcar_complete(directory):
                return JobStatus(JobState.COMPLETE)
            else:
                error = self._check_outcar_error(directory)
                return JobStatus(
                    JobState.FAILED,
                    message=error or f"Exit code {result.returncode}"
                )
        except Exception as e:
            return JobStatus(JobState.FAILED, message=str(e))

    def _run_background(self, directory: str, cmd: str) -> int:
        """Start VASP in background and return PID."""
        pid_file = os.path.join(directory, '.vasp_pid')
        log_file = os.path.join(directory, 'vasp.log')

        # Use nohup to survive terminal close
        full_cmd = f"cd {directory} && nohup {cmd} > {log_file} 2>&1 & echo $!"

        result = subprocess.run(
            full_cmd,
            shell=True,
            capture_output=True,
            text=True,
        )

        pid = int(result.stdout.strip())
        with open(pid_file, 'w') as f:
            f.write(str(pid))

        return pid

    def _read_pid(self, directory: str) -> int | None:
        """Read PID from tracking file."""
        pid_file = os.path.join(directory, '.vasp_pid')
        if os.path.exists(pid_file):
            try:
                with open(pid_file) as f:
                    return int(f.read().strip())
            except (OSError, ValueError):
                return None
        return None

    def _is_process_running(self, pid: int) -> bool:
        """Check if process with given PID is running."""
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False

    def _verify_inputs(self, directory: str) -> None:
        """Verify required input files exist."""
        required = ['INCAR', 'POSCAR', 'POTCAR', 'KPOINTS']
        missing = []

        for fname in required:
            if not os.path.exists(os.path.join(directory, fname)):
                # KPOINTS not required if KSPACING is set
                if fname == 'KPOINTS':
                    incar = os.path.join(directory, 'INCAR')
                    if os.path.exists(incar):
                        with open(incar) as f:
                            if 'KSPACING' in f.read().upper():
                                continue
                missing.append(fname)

        if missing:
            raise VaspSetupError(
                f"Missing input files in {directory}: {', '.join(missing)}"
            )

    def __repr__(self) -> str:
        return (
            f"LocalRunner(vasp_command={self.vasp_command!r}, "
            f"mpi_command={self.mpi_command!r}, background={self.background})"
        )
