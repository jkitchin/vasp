"""Base runner interface for VASP execution.

This module defines the abstract Runner interface that all
execution backends must implement. Runners are designed for
async/non-blocking operation suitable for long-running calculations.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass


class JobState(Enum):
    """Possible states of a VASP calculation."""

    NOT_STARTED = "not_started"  # No job submitted, no results
    SUBMITTED = "submitted"  # Just submitted this run
    QUEUED = "queued"  # Waiting in queue
    RUNNING = "running"  # Currently executing
    COMPLETE = "complete"  # Finished successfully
    FAILED = "failed"  # Finished with error
    UNKNOWN = "unknown"  # Cannot determine state


@dataclass
class JobStatus:
    """Status information for a VASP calculation.

    Attributes:
        state: Current state of the job.
        jobid: Job identifier from the scheduler (if applicable).
        message: Additional status message or error description.
        metadata: Additional runner-specific data.
    """

    state: JobState
    jobid: str | None = None
    message: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def is_done(self) -> bool:
        """Check if calculation has finished (success or failure)."""
        return self.state in (JobState.COMPLETE, JobState.FAILED)

    @property
    def is_active(self) -> bool:
        """Check if calculation is queued or running."""
        return self.state in (JobState.QUEUED, JobState.RUNNING, JobState.SUBMITTED)

    @property
    def is_success(self) -> bool:
        """Check if calculation completed successfully."""
        return self.state == JobState.COMPLETE


class Runner(ABC):
    """Abstract base class for VASP execution backends.

    Runners handle submitting and monitoring VASP calculations.
    They are designed for non-blocking operation:
    - run() submits/starts a job and returns immediately
    - status() checks current state without blocking
    - Exceptions signal state transitions to calling code

    Subclasses must implement:
    - run(): Submit or start a calculation
    - status(): Check current job status

    Example:
        >>> runner = LocalRunner(vasp_command='vasp_std')
        >>> status = runner.run('/path/to/calc')
        >>> if status.state == JobState.COMPLETE:
        ...     print("Done!")
    """

    @abstractmethod
    def run(self, directory: str) -> JobStatus:
        """Start or submit a VASP calculation.

        This method should NOT block waiting for completion.
        For async runners (SLURM, K8s), it submits the job.
        For local runners, behavior depends on configuration.

        Args:
            directory: Path to calculation directory containing input files.

        Returns:
            JobStatus with current state after submission.

        Raises:
            VaspSubmitted: Job was just submitted (contains jobid).
            VaspQueued: Job is already queued.
            VaspRunning: Job is currently running.
            VaspSetupError: Input files missing or invalid.
        """
        pass

    @abstractmethod
    def status(self, directory: str) -> JobStatus:
        """Check status of a calculation without blocking.

        This method never starts a new calculation. It only
        reports the current state based on job scheduler
        status and output files.

        Args:
            directory: Path to calculation directory.

        Returns:
            JobStatus with current state.
        """
        pass

    def cancel(self, directory: str) -> bool:
        """Cancel a running or queued job.

        Args:
            directory: Path to calculation directory.

        Returns:
            True if cancellation was successful or job wasn't running.
        """
        return False

    def wait(
        self, directory: str, timeout: float | None = None, poll_interval: float = 30.0
    ) -> JobStatus:
        """Block until calculation completes (for interactive use).

        This is a convenience method that polls status() until
        the job is done. For most workflows, you should use
        status() directly with your own retry logic.

        Args:
            directory: Path to calculation directory.
            timeout: Maximum seconds to wait (None = forever).
            poll_interval: Seconds between status checks.

        Returns:
            Final JobStatus.

        Raises:
            TimeoutError: If timeout is reached before completion.
        """
        import time

        start = time.time()
        while True:
            s = self.status(directory)
            if s.is_done:
                return s
            if timeout is not None and (time.time() - start) > timeout:
                raise TimeoutError(f"Calculation not complete after {timeout}s")
            time.sleep(poll_interval)

    def get_logs(self, directory: str, tail_lines: int = 100) -> str:
        """Get log output from the calculation.

        Args:
            directory: Path to calculation directory.
            tail_lines: Number of lines to return from end of log.

        Returns:
            Log content as string.
        """
        import os

        # Default implementation reads OUTCAR
        outcar = os.path.join(directory, "OUTCAR")
        if os.path.exists(outcar):
            with open(outcar) as f:
                lines = f.readlines()
                return "".join(lines[-tail_lines:])
        return "No OUTCAR found"

    def _check_outcar_complete(self, directory: str) -> bool:
        """Check if OUTCAR indicates successful completion.

        Args:
            directory: Path to calculation directory.

        Returns:
            True if OUTCAR contains completion marker.
        """
        import os

        # Check main directory first
        outcar = os.path.join(directory, "OUTCAR")
        if os.path.exists(outcar):
            if self._outcar_has_completion_marker(outcar):
                return True

        # Check for NEB calculation (subdirectories 01/, 02/, etc.)
        # For NEB, check if first intermediate image (01/) has completed
        neb_outcar = os.path.join(directory, "01", "OUTCAR")
        if os.path.exists(neb_outcar):
            return self._outcar_has_completion_marker(neb_outcar)

        return False

    def _outcar_has_completion_marker(self, outcar: str) -> bool:
        """Check if an OUTCAR file has the completion marker."""
        # Check last few KB for completion marker
        with open(outcar, "rb") as f:
            f.seek(0, 2)  # End of file
            size = f.tell()
            # Read last 10KB
            f.seek(max(0, size - 10240))
            content = f.read().decode("utf-8", errors="ignore")

        return "General timing and accounting" in content

    def _check_outcar_error(self, directory: str) -> str | None:
        """Check OUTCAR for error messages.

        Args:
            directory: Path to calculation directory.

        Returns:
            Error message if found, None otherwise.
        """
        import os

        outcar = os.path.join(directory, "OUTCAR")
        if not os.path.exists(outcar):
            return None

        error_patterns = [
            "ZBRENT: fatal error",
            "VERY BAD NEWS!",
            "internal error",
            "EDDDAV: Call to ZHEGV failed",
            "Sub-Space-Matrix is not hermitian",
        ]

        with open(outcar, "rb") as f:
            content = f.read().decode("utf-8", errors="ignore")

        for pattern in error_patterns:
            if pattern in content:
                return pattern

        return None

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
