"""Tests for VASP exceptions."""

import pytest

from vasp.exceptions import (
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


class TestVaspException:
    """Test base VaspException."""

    def test_basic_exception(self):
        """Test basic exception creation."""
        exc = VaspException("Something went wrong")

        assert str(exc) == "Something went wrong"

    def test_exception_inheritance(self):
        """Test that VaspException inherits from Exception."""
        exc = VaspException("test")

        assert isinstance(exc, Exception)


class TestVaspSubmitted:
    """Test VaspSubmitted exception."""

    def test_with_jobid(self):
        """Test exception with job ID."""
        exc = VaspSubmitted(jobid="12345")

        assert exc.jobid == "12345"
        assert "12345" in str(exc)

    def test_with_message(self):
        """Test exception with custom message."""
        exc = VaspSubmitted(message="Job submitted to SLURM", jobid="99999")

        assert "submitted" in str(exc).lower()
        assert exc.jobid == "99999"

    def test_raise_and_catch(self):
        """Test raising and catching the exception."""
        with pytest.raises(VaspSubmitted) as exc_info:
            raise VaspSubmitted(jobid="abc123")

        assert exc_info.value.jobid == "abc123"


class TestVaspQueued:
    """Test VaspQueued exception."""

    def test_with_jobid(self):
        """Test exception with job ID."""
        exc = VaspQueued(jobid="54321")

        assert exc.jobid == "54321"

    def test_catch_as_base(self):
        """Test catching as VaspException."""
        with pytest.raises(VaspException):
            raise VaspQueued(jobid="test")


class TestVaspRunning:
    """Test VaspRunning exception."""

    def test_with_message(self):
        """Test exception with progress message."""
        exc = VaspRunning(message="Iteration 5 of 100", jobid="run-1")

        assert exc.jobid == "run-1"
        assert "Iteration" in str(exc)

    def test_default_message(self):
        """Test default message."""
        exc = VaspRunning(jobid="xyz")

        assert "running" in str(exc).lower()


class TestVaspNotFinished:
    """Test VaspNotFinished exception."""

    def test_basic(self):
        """Test basic exception."""
        exc = VaspNotFinished("Calculation incomplete")

        assert "incomplete" in str(exc).lower()


class TestVaspNotConverged:
    """Test VaspNotConverged exception."""

    def test_with_details(self):
        """Test with convergence details."""
        exc = VaspNotConverged("Electronic SCF did not converge after 60 steps")

        assert "converge" in str(exc).lower()
        assert "60" in str(exc)


class TestVaspError:
    """Test VaspError exception."""

    def test_basic(self):
        """Test basic error."""
        exc = VaspError("POTCAR mismatch")

        assert "POTCAR" in str(exc)


class TestVaspEmptyOutput:
    """Test VaspEmptyOutput exception."""

    def test_basic(self):
        """Test empty output error."""
        exc = VaspEmptyOutput("OUTCAR is empty")

        assert "empty" in str(exc).lower()


class TestVaspSetupError:
    """Test VaspSetupError exception."""

    def test_missing_files(self):
        """Test missing files error."""
        exc = VaspSetupError("Missing POTCAR for element Fe")

        assert "POTCAR" in str(exc)
        assert "Fe" in str(exc)


class TestVaspWarning:
    """Test VaspWarning."""

    def test_basic(self):
        """Test basic warning."""
        with pytest.warns(VaspWarning):
            import warnings
            warnings.warn("Check ENCUT setting", VaspWarning, stacklevel=2)


class TestExceptionHierarchy:
    """Test exception hierarchy and catching."""

    def test_catch_submitted_as_exception(self):
        """Test catching submitted as VaspException."""
        with pytest.raises(VaspException):
            raise VaspSubmitted(jobid="test")

    def test_catch_queued_as_exception(self):
        """Test catching queued as VaspException."""
        with pytest.raises(VaspException):
            raise VaspQueued(jobid="test")

    def test_catch_running_as_exception(self):
        """Test catching running as VaspException."""
        with pytest.raises(VaspException):
            raise VaspRunning(jobid="test")

    def test_distinguish_exceptions(self):
        """Test distinguishing between exception types."""
        exceptions = [
            VaspSubmitted(jobid="1"),
            VaspQueued(jobid="2"),
            VaspRunning(jobid="3"),
        ]

        submitted_count = sum(1 for e in exceptions if isinstance(e, VaspSubmitted))
        queued_count = sum(1 for e in exceptions if isinstance(e, VaspQueued))
        running_count = sum(1 for e in exceptions if isinstance(e, VaspRunning))

        assert submitted_count == 1
        assert queued_count == 1
        assert running_count == 1

    def test_exception_in_workflow(self):
        """Test exception handling in typical workflow pattern."""
        states = []

        for _ in range(4):
            try:
                # Simulate different states
                if len(states) == 0:
                    raise VaspSubmitted(jobid="job-1")
                elif len(states) == 1:
                    raise VaspQueued(jobid="job-1")
                elif len(states) == 2:
                    raise VaspRunning(jobid="job-1")
                else:
                    states.append("complete")
                    break
            except VaspSubmitted:
                states.append("submitted")
            except VaspQueued:
                states.append("queued")
            except VaspRunning:
                states.append("running")

        assert states == ["submitted", "queued", "running", "complete"]
