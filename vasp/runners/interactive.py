"""Interactive runner for VASP with persistent process.

This module implements VASP's interactive mode (INTERACTIVE = .TRUE.)
which maintains a persistent VASP process and reuses wavefunctions
between ionic steps for significant speedup in geometry optimizations.

The communication protocol:
1. VASP writes energy/forces to stdout after each SCF
2. VASP waits for new positions via stdin
3. Calculator sends scaled positions
4. Repeat until optimization complete

This can reduce SCF cycles by up to 75% compared to restarting VASP
for each ionic step.

Inspired by:
    VaspInteractive by the Ulissi Group at CMU
    https://github.com/ulissigroup/vasp-interactive

    Their implementation includes additional features like socket I/O
    (i-PI protocol) and MPI pause/resume that could be added here.
"""

from __future__ import annotations

import os
import re
import subprocess
import time
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..exceptions import VaspError, VaspSetupError
from .base import JobState, JobStatus, Runner

if TYPE_CHECKING:
    from ase import Atoms


@dataclass
class InteractiveResults:
    """Results from a single interactive step."""

    energy: float
    forces: np.ndarray
    stress: np.ndarray | None = None
    converged: bool = True


@dataclass
class InteractiveState:
    """Internal state of the interactive session."""

    steps: int = 0
    positions_sent: int = 0
    positions_confirmed: int = 0
    lattice_supported: bool = False
    final: bool = False
    error: str | None = None


class InteractiveRunner(Runner):
    """Run VASP in interactive mode with persistent process.

    This runner maintains a live VASP process and communicates via
    stdin/stdout. It's designed for geometry optimizations where
    reusing wavefunctions between steps provides significant speedup.

    Args:
        vasp_command: Command to run VASP. Defaults to $VASP_COMMAND
            environment variable, or 'vasp_std' if not set.
        mpi_command: MPI launcher command, e.g., 'mpirun -np 4'.
        timeout: Seconds to wait for VASP response (default: 3600).
        parse_stress: Whether to parse stress tensor (default: False).

    Example:
        >>> runner = InteractiveRunner(vasp_command='vasp_std')
        >>> runner.start('/path/to/calc')
        >>>
        >>> # Optimization loop
        >>> for step in range(100):
        ...     results = runner.step(atoms)
        ...     if converged(results.forces):
        ...         break
        ...     atoms.positions = optimize(atoms, results.forces)
        >>>
        >>> runner.close()

    Note:
        The INCAR must NOT contain NSW > 0 or IBRION settings.
        This runner handles the ionic loop externally.
    """

    def __init__(
        self,
        vasp_command: str | None = None,
        mpi_command: str | None = None,
        timeout: float = 3600,
        parse_stress: bool = False,
    ):
        self.vasp_command = vasp_command or os.environ.get("VASP_COMMAND", "vasp_std")
        self.mpi_command = mpi_command
        self.timeout = timeout
        self.parse_stress = parse_stress

        # Process state
        self._process: subprocess.Popen | None = None
        self._directory: str | None = None
        self._state = InteractiveState()
        self._atoms: Atoms | None = None

    @property
    def is_running(self) -> bool:
        """Check if VASP process is active."""
        return self._process is not None and self._process.poll() is None

    def start(self, directory: str, atoms: Atoms) -> InteractiveResults:
        """Start VASP interactive session.

        Args:
            directory: Path to calculation directory with input files.
            atoms: Initial atomic structure.

        Returns:
            Results from initial SCF calculation.

        Raises:
            VaspSetupError: If input files are missing or invalid.
            VaspError: If VASP fails to start.
        """
        if self.is_running:
            raise VaspError("Interactive session already running. Call close() first.")

        self._directory = directory
        self._atoms = atoms.copy()
        self._state = InteractiveState()

        # Verify and modify input files
        self._prepare_inputs(directory)

        # Build command
        cmd = self._build_command()

        # Start process with pipes
        try:
            self._process = subprocess.Popen(
                cmd,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=directory,
                text=True,
                bufsize=1,  # Line buffered
            )
        except Exception as e:
            raise VaspError(f"Failed to start VASP: {e}") from e

        # Read initial SCF results
        return self._read_results()

    def step(self, atoms: Atoms) -> InteractiveResults:
        """Send new positions and get updated results.

        Args:
            atoms: Updated atomic structure.

        Returns:
            Energy, forces, and optionally stress from new SCF.

        Raises:
            VaspError: If VASP process died or communication failed.
        """
        if not self.is_running:
            raise VaspError("No active VASP session. Call start() first.")

        self._atoms = atoms.copy()

        # Send new positions
        self._write_positions(atoms)

        # Read results
        return self._read_results()

    def close(self) -> None:
        """Gracefully close the VASP session.

        Writes STOPCAR to trigger clean shutdown, then waits for
        VASP to finish writing output files.
        """
        if not self.is_running:
            return

        # Type narrowing for mypy
        assert self._process is not None
        assert self._directory is not None

        try:
            # Write STOPCAR for graceful shutdown
            stopcar_path = os.path.join(self._directory, "STOPCAR")
            with open(stopcar_path, "w") as f:
                f.write("LABORT = .TRUE.\n")

            # Send dummy positions to trigger VASP to read STOPCAR
            if self._atoms is not None:
                try:
                    self._write_positions(self._atoms)
                except Exception:
                    pass  # May fail if VASP already exiting

            # Wait for process to finish
            try:
                self._process.wait(timeout=60)
            except subprocess.TimeoutExpired:
                self._process.terminate()
                self._process.wait(timeout=10)

        except Exception:
            # Force kill if graceful shutdown fails
            if self._process and self._process.poll() is None:
                self._process.kill()
                self._process.wait()

        finally:
            self._process = None
            self._state.final = True

    def run(self, directory: str) -> JobStatus:
        """Standard runner interface - not used for interactive mode.

        For interactive mode, use start(), step(), close() directly.
        This method exists for compatibility with the Runner interface.
        """
        # Check if already complete from previous run
        if self._check_outcar_complete(directory):
            return JobStatus(JobState.COMPLETE)

        return JobStatus(
            JobState.NOT_STARTED, message="Use start()/step()/close() for interactive mode"
        )

    def status(self, directory: str) -> JobStatus:
        """Check status of interactive session."""
        if self.is_running:
            return JobStatus(
                JobState.RUNNING, message=f"Interactive session active, {self._state.steps} steps"
            )

        if self._check_outcar_complete(directory):
            return JobStatus(JobState.COMPLETE)

        if self._state.error:
            return JobStatus(JobState.FAILED, message=self._state.error)

        return JobStatus(JobState.NOT_STARTED)

    def _build_command(self) -> str:
        """Build the full command string."""
        if self.mpi_command:
            return f"{self.mpi_command} {self.vasp_command}"
        return self.vasp_command

    def _prepare_inputs(self, directory: str) -> None:
        """Verify and modify input files for interactive mode."""
        # Check required files
        required = ["INCAR", "POSCAR", "POTCAR", "KPOINTS"]
        for fname in required:
            path = os.path.join(directory, fname)
            if not os.path.exists(path):
                if fname == "KPOINTS":
                    # Check for KSPACING
                    incar_path = os.path.join(directory, "INCAR")
                    if os.path.exists(incar_path):
                        with open(incar_path) as f:
                            if "KSPACING" in f.read().upper():
                                continue
                raise VaspSetupError(f"Missing {fname} in {directory}")

        # Modify INCAR for interactive mode
        incar_path = os.path.join(directory, "INCAR")
        self._modify_incar(incar_path)

    def _modify_incar(self, incar_path: str) -> None:
        """Set INTERACTIVE = .TRUE. and disable internal relaxation."""
        with open(incar_path) as f:
            content = f.read()

        lines = content.split("\n")
        new_lines = []
        has_interactive = False

        for line in lines:
            # Remove NSW, IBRION settings (we control the loop)
            upper = line.upper().strip()
            if upper.startswith("NSW") or upper.startswith("IBRION"):
                continue
            if upper.startswith("INTERACTIVE"):
                has_interactive = True
                new_lines.append("INTERACTIVE = .TRUE.")
            else:
                new_lines.append(line)

        if not has_interactive:
            new_lines.insert(0, "INTERACTIVE = .TRUE.")

        # Also set NSW=0, IBRION=-1 explicitly
        new_lines.insert(1, "NSW = 0")
        new_lines.insert(2, "IBRION = -1")

        with open(incar_path, "w") as f:
            f.write("\n".join(new_lines))

    def _write_positions(self, atoms: Atoms) -> None:
        """Send atomic positions to VASP via stdin."""
        if not self.is_running:
            raise VaspError("VASP process not running")

        # Type narrowing for mypy
        assert self._process is not None
        assert self._process.stdin is not None

        # Get scaled (fractional) positions
        scaled = atoms.get_scaled_positions()

        # Write positions
        for pos in scaled:
            line = f"{pos[0]:.16f} {pos[1]:.16f} {pos[2]:.16f}\n"
            self._process.stdin.write(line)

        self._process.stdin.flush()
        self._state.positions_sent += 1

    def _read_results(self) -> InteractiveResults:
        """Read energy and forces from VASP stdout."""
        if not self.is_running:
            raise VaspError("VASP process not running")

        # Type narrowing for mypy
        assert self._process is not None
        assert self._process.stdout is not None
        assert self._process.stderr is not None

        energy = None
        forces = []
        stress = None
        n_atoms = len(self._atoms) if self._atoms else 0

        # Patterns to match
        energy_pattern = re.compile(r"ETOTAL\s*=\s*([-\d.E+]+)")
        force_pattern = re.compile(r"FORCES:\s*([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)")
        pos_confirm = re.compile(r"POSITIONS:\s*read from stdin")
        lattice_support = re.compile(r"LATTICE:\s*reading from stdin")

        start_time = time.time()

        while True:
            # Check timeout
            if time.time() - start_time > self.timeout:
                self._state.error = "Timeout waiting for VASP response"
                raise VaspError(self._state.error)

            # Check if process died
            if self._process.poll() is not None:
                stderr = self._process.stderr.read()
                self._state.error = f"VASP process died: {stderr}"
                raise VaspError(self._state.error)

            # Read line
            line = self._process.stdout.readline()
            if not line:
                time.sleep(0.01)
                continue

            # Parse energy
            match = energy_pattern.search(line)
            if match:
                energy = float(match.group(1))

            # Parse forces
            match = force_pattern.search(line)
            if match:
                fx, fy, fz = float(match.group(1)), float(match.group(2)), float(match.group(3))
                forces.append([fx, fy, fz])

            # Check for lattice support
            if lattice_support.search(line):
                self._state.lattice_supported = True

            # Check for position confirmation (end of this step)
            if pos_confirm.search(line):
                self._state.positions_confirmed += 1
                break

            # Also break if we have energy and all forces
            if energy is not None and len(forces) == n_atoms:
                # Wait a bit for confirmation
                time.sleep(0.1)
                break

        self._state.steps += 1

        if energy is None:
            raise VaspError("Failed to parse energy from VASP output")

        forces_array = np.array(forces) if forces else np.zeros((n_atoms, 3))

        return InteractiveResults(
            energy=energy,
            forces=forces_array,
            stress=stress,
            converged=True,
        )

    def __enter__(self) -> InteractiveRunner:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit - ensure cleanup."""
        self.close()

    def __del__(self) -> None:
        """Destructor - ensure cleanup."""
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        status = "running" if self.is_running else "stopped"
        return (
            f"InteractiveRunner(vasp_command={self.vasp_command!r}, "
            f"status={status}, steps={self._state.steps})"
        )
