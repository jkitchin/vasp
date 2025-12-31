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


def get_optimal_nprocs(directory: str | None = None) -> int:
    """Determine optimal number of MPI processes.

    Uses heuristics based on available CPUs and system size:
    - Reads NCORE/NPAR from INCAR if available
    - Checks number of atoms from POSCAR
    - Returns reasonable default based on CPU count

    Args:
        directory: Calculation directory (optional, for reading INCAR/POSCAR)

    Returns:
        Recommended number of MPI processes
    """
    cpu_count = os.cpu_count() or 1

    # Start with half of available CPUs (leave room for system)
    # but at least 1 and cap at reasonable maximum
    nprocs = max(1, min(cpu_count // 2, 32))

    # Try to read INCAR for hints
    if directory:
        incar_path = os.path.join(directory, "INCAR")
        poscar_path = os.path.join(directory, "POSCAR")

        # Check number of atoms
        natoms = 1
        if os.path.exists(poscar_path):
            try:
                with open(poscar_path) as f:
                    lines = f.readlines()
                    # Line 6 or 7 has atom counts
                    for i in [5, 6]:
                        if i < len(lines):
                            parts = lines[i].split()
                            if parts and parts[0].isdigit():
                                natoms = sum(int(x) for x in parts)
                                break
            except (OSError, ValueError):
                pass

        # Small systems don't benefit from many processes
        # Rule of thumb: ~10-50 atoms per process is efficient
        max_efficient = max(1, natoms // 10)
        nprocs = min(nprocs, max_efficient)

        # Check NCORE setting in INCAR
        if os.path.exists(incar_path):
            try:
                with open(incar_path) as f:
                    content = f.read().upper()
                    for line in content.split("\n"):
                        if "NCORE" in line and "=" in line:
                            # NCORE should divide nprocs evenly
                            ncore = int(line.split("=")[1].split()[0])
                            # Adjust nprocs to be divisible by NCORE
                            if ncore > 1:
                                nprocs = (nprocs // ncore) * ncore
                                nprocs = max(ncore, nprocs)
                            break
            except (OSError, ValueError):
                pass

    return max(1, nprocs)


class LocalRunner(Runner):
    """Run VASP on the local machine.

    By default runs synchronously (blocking until complete).
    Use background=True for non-blocking execution.

    The runner intelligently determines the number of MPI processes:
    - If nprocs='auto' (default): detects based on CPUs and system size
    - If nprocs is an int: uses that exact number
    - Respects NCORE settings in INCAR

    Environment variables (in order of precedence):
    - VASP_EXECUTABLE: Path to VASP binary (e.g., /opt/vasp/bin/vasp_std)
    - VASP_COMMAND: Full command (legacy, e.g., 'mpirun -np 8 vasp_std')
    - VASP_NPROCS: Number of MPI processes (or 'auto')
    - VASP_MPI_EXTRA_ARGS: Extra MPI arguments (e.g., '--map-by hwthread')

    Args:
        vasp_executable: Path to VASP binary. Defaults to $VASP_EXECUTABLE,
            then extracts from $VASP_COMMAND, or 'vasp_std'.
        nprocs: Number of MPI processes. 'auto' (default) detects optimal
            count, or specify an integer. Set to 1 for serial execution.
        mpi_command: MPI launcher (default: 'mpirun'). Set to None to
            disable MPI and run serial VASP.
        mpi_extra_args: Extra arguments for MPI launcher.
        background: If True, run in background and return immediately.

    Example:
        >>> # Auto-detect optimal parallelization
        >>> runner = LocalRunner()
        >>> calc = Vasp('my_calc', runner=runner, ...)

        >>> # Specify exact process count
        >>> runner = LocalRunner(nprocs=16)

        >>> # Serial execution (no MPI)
        >>> runner = LocalRunner(nprocs=1, mpi_command=None)

        >>> # Custom MPI settings
        >>> runner = LocalRunner(
        ...     nprocs=8,
        ...     mpi_extra_args='--bind-to core --map-by socket'
        ... )
    """

    def __init__(
        self,
        vasp_executable: str | None = None,
        nprocs: int | str = "auto",
        mpi_command: str | None = "mpirun",
        mpi_extra_args: str | None = None,
        background: bool = False,
        # Legacy parameter
        vasp_command: str | None = None,
    ):
        # Handle legacy VASP_COMMAND that includes MPI
        legacy_cmd = vasp_command or os.environ.get("VASP_COMMAND", "")
        if legacy_cmd and ("mpirun" in legacy_cmd or "srun" in legacy_cmd):
            # Parse legacy format: "mpirun -np 8 ... vasp_std"
            self._use_legacy_command = True
            self._legacy_command = legacy_cmd
            # Extract executable for display
            parts = legacy_cmd.split()
            self.vasp_executable = parts[-1] if parts else "vasp_std"
            self.nprocs = nprocs
            self.mpi_command = None
            self.mpi_extra_args = None
        else:
            self._use_legacy_command = False
            self._legacy_command = None

            # Determine VASP executable
            self.vasp_executable = (
                vasp_executable or os.environ.get("VASP_EXECUTABLE") or legacy_cmd or "vasp_std"
            )

            # Process count
            env_nprocs = os.environ.get("VASP_NPROCS", "")
            if nprocs == "auto" and env_nprocs:
                if env_nprocs.lower() == "auto":
                    self.nprocs = "auto"
                else:
                    try:
                        self.nprocs = int(env_nprocs)
                    except ValueError:
                        self.nprocs = "auto"
            else:
                self.nprocs = nprocs

            # MPI settings
            self.mpi_command = mpi_command
            self.mpi_extra_args = mpi_extra_args or os.environ.get("VASP_MPI_EXTRA_ARGS")

        self.background = background
        self._current_directory: str | None = None

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

        # If previously failed, clean up old output files so we can re-run
        if current.state == JobState.FAILED:
            self._clean_output_files(directory)

        # Verify input files exist
        self._verify_inputs(directory)

        # Build command (pass directory for auto nprocs detection)
        cmd = self._build_command(directory)

        if self.background:
            pid = self._run_background(directory, cmd)
            raise VaspRunning(message=f"Started background process {pid}", jobid=str(pid))
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
        outcar = os.path.join(directory, "OUTCAR")
        if os.path.exists(outcar):
            # OUTCAR exists but not complete and no process
            return JobStatus(JobState.FAILED, message="OUTCAR incomplete and no running process")

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

    def _build_command(self, directory: str) -> str:
        """Build the full command string.

        Args:
            directory: Calculation directory (used for auto nprocs detection)

        Returns:
            Full command string to execute
        """
        # Use legacy command if set
        if self._use_legacy_command and self._legacy_command:
            return self._legacy_command

        # Determine number of processes
        if self.nprocs == "auto":
            nprocs = get_optimal_nprocs(directory)
        else:
            nprocs = int(self.nprocs)

        # Build MPI command if needed
        if self.mpi_command and nprocs > 1:
            parts = [self.mpi_command, "-np", str(nprocs)]
            if self.mpi_extra_args:
                parts.append(self.mpi_extra_args)
            parts.append(self.vasp_executable)
            return " ".join(parts)

        # Serial execution
        return self.vasp_executable

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
                return JobStatus(JobState.FAILED, message=error or f"Exit code {result.returncode}")
        except Exception as e:
            return JobStatus(JobState.FAILED, message=str(e))

    def _run_background(self, directory: str, cmd: str) -> int:
        """Start VASP in background and return PID."""
        pid_file = os.path.join(directory, ".vasp_pid")
        log_file = os.path.join(directory, "vasp.log")

        # Use nohup to survive terminal close
        full_cmd = f"cd {directory} && nohup {cmd} > {log_file} 2>&1 & echo $!"

        result = subprocess.run(
            full_cmd,
            shell=True,
            capture_output=True,
            text=True,
        )

        pid = int(result.stdout.strip())
        with open(pid_file, "w") as f:
            f.write(str(pid))

        return pid

    def _read_pid(self, directory: str) -> int | None:
        """Read PID from tracking file."""
        pid_file = os.path.join(directory, ".vasp_pid")
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

    def _clean_output_files(self, directory: str) -> None:
        """Remove output files from a failed run to allow re-running."""
        output_files = [
            "OUTCAR",
            "OSZICAR",
            "vasprun.xml",
            "CONTCAR",
            "CHG",
            "CHGCAR",
            "WAVECAR",
            "PROCAR",
            "EIGENVAL",
            "DOSCAR",
            "PCDAT",
            "XDATCAR",
            "REPORT",
            ".vasp_pid",
        ]
        for fname in output_files:
            fpath = os.path.join(directory, fname)
            if os.path.exists(fpath):
                os.remove(fpath)

    def _verify_inputs(self, directory: str) -> None:
        """Verify required input files exist."""
        required = ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]
        missing = []

        # Check if this is a NEB calculation (IMAGES in INCAR or 00/ subdir exists)
        is_neb = False
        incar = os.path.join(directory, "INCAR")
        if os.path.exists(incar):
            with open(incar) as f:
                if "IMAGES" in f.read().upper():
                    is_neb = True
        if os.path.exists(os.path.join(directory, "00", "POSCAR")):
            is_neb = True

        for fname in required:
            if not os.path.exists(os.path.join(directory, fname)):
                # POSCAR not required in main dir for NEB (it's in subdirs)
                if fname == "POSCAR" and is_neb:
                    continue
                # KPOINTS not required if KSPACING is set
                if fname == "KPOINTS":
                    if os.path.exists(incar):
                        with open(incar) as f:
                            if "KSPACING" in f.read().upper():
                                continue
                missing.append(fname)

        if missing:
            raise VaspSetupError(f"Missing input files in {directory}: {', '.join(missing)}")

    def __repr__(self) -> str:
        if self._use_legacy_command:
            return f"LocalRunner(vasp_command={self._legacy_command!r})"
        return (
            f"LocalRunner(vasp_executable={self.vasp_executable!r}, "
            f"nprocs={self.nprocs!r}, mpi_command={self.mpi_command!r})"
        )
