"""Tests for InteractiveRunner."""

import os
from unittest.mock import MagicMock

import numpy as np
import pytest

from vasp.exceptions import VaspError, VaspSetupError
from vasp.runners import InteractiveResults, InteractiveRunner
from vasp.runners.interactive import InteractiveState

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def temp_calc_dir(tmp_path):
    """Create a temporary calculation directory with input files."""
    calc_dir = tmp_path / "calc"
    calc_dir.mkdir()

    # Create minimal input files
    (calc_dir / "INCAR").write_text("ENCUT = 400\nISMEAR = 0\n")
    (calc_dir / "POSCAR").write_text(
        """Si
1.0
5.43 0.0 0.0
0.0 5.43 0.0
0.0 0.0 5.43
Si
2
Direct
0.0 0.0 0.0
0.25 0.25 0.25
"""
    )
    (calc_dir / "POTCAR").write_text("PAW_PBE Si\n")
    (calc_dir / "KPOINTS").write_text("Automatic\n0\nGamma\n4 4 4\n")

    return str(calc_dir)


@pytest.fixture
def si_atoms():
    """Create Si atoms for testing."""
    from ase.build import bulk

    return bulk("Si", "diamond", a=5.43)


@pytest.fixture
def mock_vasp_output():
    """Mock VASP interactive output."""
    return """
 running on    4 total cores
 distrk:  each k-point on    4 cores,    1 groups
 vasp.6.3.0 13Jan22 (build Jan 13 2022 12:52:14) complex

 POSCAR found :  2 types and       2 ions
 LDA part: xc-table for Pade apance of LDA

 POTCAR:    PAW_PBE Si 05Jan2001

 ETOTAL = -10.84274516
 FORCES:  0.00000000  0.00000000  0.00000000
 FORCES:  0.00000000  0.00000000  0.00000000
 POSITIONS: read from stdin
"""


# =============================================================================
# InteractiveResults Tests
# =============================================================================


class TestInteractiveResults:
    """Tests for InteractiveResults dataclass."""

    def test_create_results(self):
        """Test creating InteractiveResults."""
        forces = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]])
        results = InteractiveResults(
            energy=-10.5,
            forces=forces,
            stress=None,
            converged=True,
        )

        assert results.energy == -10.5
        assert results.forces.shape == (2, 3)
        assert results.converged is True

    def test_results_with_stress(self):
        """Test results with stress tensor."""
        stress = np.array([1.0, 2.0, 3.0, 0.1, 0.2, 0.3])
        results = InteractiveResults(
            energy=-5.0,
            forces=np.zeros((4, 3)),
            stress=stress,
        )

        assert results.stress is not None
        assert len(results.stress) == 6


# =============================================================================
# InteractiveState Tests
# =============================================================================


class TestInteractiveState:
    """Tests for InteractiveState dataclass."""

    def test_initial_state(self):
        """Test initial state values."""
        state = InteractiveState()

        assert state.steps == 0
        assert state.positions_sent == 0
        assert state.positions_confirmed == 0
        assert state.lattice_supported is False
        assert state.final is False
        assert state.error is None


# =============================================================================
# InteractiveRunner Tests
# =============================================================================


class TestInteractiveRunner:
    """Tests for InteractiveRunner class."""

    def test_init_defaults(self, monkeypatch):
        """Test default initialization."""
        # Clear VASP_COMMAND env var for this test
        monkeypatch.delenv("VASP_COMMAND", raising=False)

        runner = InteractiveRunner()

        assert runner.vasp_command == "vasp_std"
        assert runner.mpi_command is None
        assert runner.timeout == 3600
        assert runner.parse_stress is False
        assert not runner.is_running

    def test_init_custom(self):
        """Test custom initialization."""
        runner = InteractiveRunner(
            vasp_command="vasp_gam",
            mpi_command="mpirun -np 8",
            timeout=7200,
            parse_stress=True,
        )

        assert runner.vasp_command == "vasp_gam"
        assert runner.mpi_command == "mpirun -np 8"
        assert runner.timeout == 7200
        assert runner.parse_stress is True

    def test_build_command_no_mpi(self):
        """Test command building without MPI."""
        runner = InteractiveRunner(vasp_command="vasp_std")
        cmd = runner._build_command()

        assert cmd == "vasp_std"

    def test_build_command_with_mpi(self):
        """Test command building with MPI."""
        runner = InteractiveRunner(
            vasp_command="vasp_std",
            mpi_command="mpirun -np 4",
        )
        cmd = runner._build_command()

        assert cmd == "mpirun -np 4 vasp_std"

    def test_prepare_inputs_missing_files(self, tmp_path):
        """Test error on missing input files."""
        runner = InteractiveRunner()
        empty_dir = str(tmp_path / "empty")
        os.makedirs(empty_dir)

        with pytest.raises(VaspSetupError, match="Missing INCAR"):
            runner._prepare_inputs(empty_dir)

    def test_prepare_inputs_modifies_incar(self, temp_calc_dir):
        """Test that INCAR is modified for interactive mode."""
        runner = InteractiveRunner()
        runner._prepare_inputs(temp_calc_dir)

        with open(os.path.join(temp_calc_dir, "INCAR")) as f:
            content = f.read()

        assert "INTERACTIVE = .TRUE." in content
        assert "NSW = 0" in content
        assert "IBRION = -1" in content

    def test_prepare_inputs_removes_nsw_ibrion(self, temp_calc_dir):
        """Test that NSW and IBRION are removed from INCAR."""
        # Add NSW and IBRION to INCAR
        incar_path = os.path.join(temp_calc_dir, "INCAR")
        with open(incar_path, "w") as f:
            f.write("ENCUT = 400\nNSW = 100\nIBRION = 2\nISMEAR = 0\n")

        runner = InteractiveRunner()
        runner._prepare_inputs(temp_calc_dir)

        with open(incar_path) as f:
            content = f.read()

        # Should have removed user NSW/IBRION and set our own
        assert content.count("NSW") == 1  # Only our NSW = 0
        assert content.count("IBRION") == 1  # Only our IBRION = -1

    def test_repr_stopped(self):
        """Test string representation when stopped."""
        runner = InteractiveRunner()
        repr_str = repr(runner)

        assert "InteractiveRunner" in repr_str
        assert "vasp_std" in repr_str
        assert "stopped" in repr_str

    def test_is_running_no_process(self):
        """Test is_running when no process exists."""
        runner = InteractiveRunner()
        assert not runner.is_running

    def test_step_without_start_raises(self, si_atoms):
        """Test that step() raises without start()."""
        runner = InteractiveRunner()

        with pytest.raises(VaspError, match="No active VASP session"):
            runner.step(si_atoms)

    def test_context_manager(self):
        """Test context manager protocol."""
        runner = InteractiveRunner()

        with runner as r:
            assert r is runner

        # close() should have been called (no error even without process)

    def test_status_not_started(self, temp_calc_dir):
        """Test status when not started."""
        runner = InteractiveRunner()
        status = runner.status(temp_calc_dir)

        assert status.state.value == "not_started"


class TestInteractiveRunnerMocked:
    """Tests with mocked subprocess for full workflow."""

    def test_write_positions(self, si_atoms):
        """Test writing positions to stdin."""
        runner = InteractiveRunner()
        runner._atoms = si_atoms

        # Create mock process
        mock_stdin = MagicMock()
        runner._process = MagicMock()
        runner._process.stdin = mock_stdin
        runner._process.poll.return_value = None  # Process running

        runner._write_positions(si_atoms)

        # Should have written 2 lines (2 atoms)
        assert mock_stdin.write.call_count == 2
        mock_stdin.flush.assert_called_once()

    def test_close_writes_stopcar(self, temp_calc_dir, si_atoms):
        """Test that close() writes STOPCAR."""
        runner = InteractiveRunner()
        runner._directory = temp_calc_dir
        runner._atoms = si_atoms

        # Create mock process
        runner._process = MagicMock()
        runner._process.poll.return_value = None  # Process running initially
        runner._process.stdin = MagicMock()
        runner._process.wait.return_value = 0

        # Make poll return 0 after first check (process finished)
        runner._process.poll.side_effect = [None, 0]

        runner.close()

        # STOPCAR should exist
        stopcar_path = os.path.join(temp_calc_dir, "STOPCAR")
        assert os.path.exists(stopcar_path)
        with open(stopcar_path) as f:
            assert "LABORT" in f.read()

    def test_start_already_running_raises(self, temp_calc_dir, si_atoms):
        """Test that start() raises if already running."""
        runner = InteractiveRunner()

        # Mock a running process
        runner._process = MagicMock()
        runner._process.poll.return_value = None

        with pytest.raises(VaspError, match="already running"):
            runner.start(temp_calc_dir, si_atoms)


class TestInteractiveResultsParsing:
    """Tests for parsing VASP output."""

    def test_parse_energy_pattern(self):
        """Test energy regex pattern."""
        import re

        pattern = re.compile(r"ETOTAL\s*=\s*([-\d.E+]+)")

        line = " ETOTAL = -10.84274516"
        match = pattern.search(line)
        assert match is not None
        assert float(match.group(1)) == pytest.approx(-10.84274516)

    def test_parse_forces_pattern(self):
        """Test forces regex pattern."""
        import re

        pattern = re.compile(r"FORCES:\s*([-\d.E+]+)\s+([-\d.E+]+)\s+([-\d.E+]+)")

        line = " FORCES:  0.12345678  -0.23456789  0.34567890"
        match = pattern.search(line)
        assert match is not None
        assert float(match.group(1)) == pytest.approx(0.12345678)
        assert float(match.group(2)) == pytest.approx(-0.23456789)
        assert float(match.group(3)) == pytest.approx(0.34567890)

    def test_parse_position_confirmation(self):
        """Test position confirmation regex."""
        import re

        pattern = re.compile(r"POSITIONS:\s*read from stdin")

        line = " POSITIONS: read from stdin"
        assert pattern.search(line) is not None

        line = " LATTICE: reading from stdin"
        assert pattern.search(line) is None
