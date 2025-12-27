"""Exceptions for VASP calculator.

These exceptions are used for non-blocking async workflow patterns.
When a calculation is not complete, an exception is raised to signal
the current state, allowing scripts to handle job submission and
monitoring without blocking.

Example:
    >>> from vasp import Vasp
    >>> from vasp.exceptions import VaspSubmitted, VaspQueued
    >>>
    >>> calc = Vasp('my_calc', atoms=atoms, runner=runner)
    >>> try:
    ...     energy = calc.potential_energy
    ... except VaspSubmitted as e:
    ...     print(f"Job submitted: {e.jobid}")
    ... except VaspQueued:
    ...     print("Job in queue, check back later")
"""

from __future__ import annotations


class VaspException(Exception):
    """Base exception for all VASP-related errors."""
    pass


class VaspSubmitted(VaspException):
    """Raised when a job has just been submitted.

    Attributes:
        jobid: The job identifier from the scheduler.
        message: Additional status message.
    """

    def __init__(self, message: str = "Job submitted", jobid: str | None = None):
        self.message = message
        self.jobid = jobid
        super().__init__(message)

    def __str__(self) -> str:
        if self.jobid:
            return f"{self.message} ({self.jobid})"
        return self.message


class VaspQueued(VaspException):
    """Raised when a job is waiting in the queue.

    Attributes:
        jobid: The job identifier.
        message: Additional status message.
    """

    def __init__(self, message: str = "Queued", jobid: str | None = None):
        self.message = message
        self.jobid = jobid
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspRunning(VaspException):
    """Raised when a job is currently executing.

    Attributes:
        jobid: The job identifier.
        message: Additional status message.
    """

    def __init__(self, message: str = "Running", jobid: str | None = None):
        self.message = message
        self.jobid = jobid
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspNotFinished(VaspException):
    """Raised when a calculation started but is not complete.

    This can indicate the job is still running or was interrupted.
    """

    def __init__(self, message: str = "Calculation not finished"):
        self.message = message
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspNotConverged(VaspException):
    """Raised when a calculation did not converge.

    This typically means the electronic or ionic relaxation
    did not reach the specified convergence criteria.
    """

    def __init__(self, message: str = "Calculation did not converge"):
        self.message = message
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspError(VaspException):
    """Raised when VASP encounters an error during calculation.

    This indicates a fatal error in the VASP run, such as
    incorrect input parameters or numerical instabilities.
    """

    def __init__(self, message: str = "VASP error"):
        self.message = message
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspEmptyOutput(VaspException):
    """Raised when expected output files are empty or missing.

    This can happen if VASP crashed before writing output
    or if the CONTCAR is empty after a failed relaxation.
    """

    def __init__(self, message: str = "Empty or missing output"):
        self.message = message
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspSetupError(VaspException):
    """Raised when there's an error in calculation setup.

    This includes missing POTCAR files, invalid parameters,
    or configuration errors.
    """

    def __init__(self, message: str = "Setup error"):
        self.message = message
        super().__init__(message)

    def __str__(self) -> str:
        return self.message


class VaspWarning(UserWarning):
    """Warning for non-fatal issues that may affect results.

    These are issues that don't prevent the calculation from
    completing but may indicate problems with the results.
    """
    pass


# Backward compatibility aliases
VaspUnknownState = VaspException
VaspEmptyCONTCAR = VaspEmptyOutput
