"""Tests for VASP runners."""

import os

import numpy as np
import pytest

from vasp.exceptions import VaspRunning, VaspSetupError, VaspSubmitted
from vasp.runners import (
    JobState,
    JobStatus,
    KubernetesRunner,
    LocalRunner,
    MockResults,
    MockRunner,
    SlurmRunner,
)


class TestJobStatus:
    """Test JobStatus dataclass."""

    def test_is_done(self):
        """Test is_done property."""
        assert JobStatus(JobState.COMPLETE).is_done
        assert JobStatus(JobState.FAILED).is_done
        assert not JobStatus(JobState.RUNNING).is_done
        assert not JobStatus(JobState.QUEUED).is_done

    def test_is_active(self):
        """Test is_active property."""
        assert JobStatus(JobState.RUNNING).is_active
        assert JobStatus(JobState.QUEUED).is_active
        assert JobStatus(JobState.SUBMITTED).is_active
        assert not JobStatus(JobState.COMPLETE).is_active

    def test_is_success(self):
        """Test is_success property."""
        assert JobStatus(JobState.COMPLETE).is_success
        assert not JobStatus(JobState.FAILED).is_success

    def test_with_jobid(self):
        """Test JobStatus with job ID."""
        status = JobStatus(JobState.RUNNING, jobid="12345")
        assert status.jobid == "12345"

    def test_with_message(self):
        """Test JobStatus with message."""
        status = JobStatus(JobState.FAILED, message="Out of memory")
        assert status.message == "Out of memory"


class TestMockRunner:
    """Test MockRunner for testing."""

    def test_basic_mock(self, temp_dir):
        """Test basic mock runner behavior."""
        results = MockResults(energy=-10.0)
        runner = MockRunner(results=results)

        status = runner.run(temp_dir)

        assert status.state == JobState.COMPLETE

    def test_mock_writes_files(self, temp_dir):
        """Test that mock runner writes output files."""
        results = MockResults(energy=-10.0, forces=np.array([[0.1, -0.1, 0.0]]), fermi_level=5.0)
        runner = MockRunner(results=results)

        runner.run(temp_dir)

        assert os.path.exists(os.path.join(temp_dir, "OUTCAR"))
        assert os.path.exists(os.path.join(temp_dir, "vasprun.xml"))

    def test_mock_state_sequence(self, temp_dir):
        """Test mock runner with state sequence."""
        runner = MockRunner(
            results=MockResults(energy=-10.0),
            state_sequence=[JobState.SUBMITTED, JobState.RUNNING, JobState.COMPLETE],
        )

        # First call - SUBMITTED (raises exception)
        with pytest.raises(VaspSubmitted):
            runner.run(temp_dir)

        # Second call - RUNNING
        status = runner.status(temp_dir)
        assert status.state == JobState.RUNNING

        # Third call - COMPLETE
        status = runner.status(temp_dir)
        assert status.state == JobState.COMPLETE

    def test_mock_delay(self, temp_dir):
        """Test mock runner with delay."""
        import time

        runner = MockRunner(results=MockResults(energy=-10.0), delay=0.1)

        start = time.time()
        runner.run(temp_dir)
        elapsed = time.time() - start

        assert elapsed >= 0.1

    def test_mock_status_not_started(self, temp_dir):
        """Test status for not started calculation."""
        runner = MockRunner(results=MockResults(energy=-10.0))

        status = runner.status(temp_dir)
        assert status.state == JobState.NOT_STARTED

    def test_mock_status_after_run(self, temp_dir):
        """Test status after running."""
        runner = MockRunner(results=MockResults(energy=-10.0))
        runner.run(temp_dir)

        status = runner.status(temp_dir)
        assert status.state == JobState.COMPLETE

    def test_mock_cancel(self, temp_dir):
        """Test cancelling mock runner."""
        runner = MockRunner(results=MockResults(energy=-10.0), state_sequence=[JobState.RUNNING])
        # run() raises VaspRunning for RUNNING state
        with pytest.raises(VaspRunning):
            runner.run(temp_dir)

        assert runner.cancel(temp_dir)


class TestGetOptimalNprocs:
    """Test get_optimal_nprocs function."""

    def test_without_directory(self):
        """Test auto-detection without directory."""
        import os

        from vasp.runners.local import get_optimal_nprocs

        nprocs = get_optimal_nprocs()
        cpu_count = os.cpu_count() or 1

        # Should be at least 1 and at most half of CPUs (capped at 32)
        assert 1 <= nprocs <= min(cpu_count // 2, 32)

    def test_with_small_system(self, temp_dir):
        """Test with small system (few atoms)."""
        from vasp.runners.local import get_optimal_nprocs

        # Create a POSCAR with 2 atoms
        poscar = os.path.join(temp_dir, "POSCAR")
        with open(poscar, "w") as f:
            f.write("Si2\n1.0\n5.43 0.0 0.0\n0.0 5.43 0.0\n0.0 0.0 5.43\n")
            f.write("Si\n2\nDirect\n0.0 0.0 0.0\n0.25 0.25 0.25\n")

        nprocs = get_optimal_nprocs(temp_dir)

        # Small system should use few processes
        assert nprocs <= 2

    def test_respects_ncore(self, temp_dir):
        """Test that NCORE in INCAR affects nprocs."""
        from vasp.runners.local import get_optimal_nprocs

        # Create POSCAR with 100 atoms
        poscar = os.path.join(temp_dir, "POSCAR")
        with open(poscar, "w") as f:
            f.write("Test\n1.0\n10 0 0\n0 10 0\n0 0 10\nSi\n100\nCartesian\n")
            for i in range(100):
                f.write(f"{i % 10} {(i // 10) % 10} {i // 100} \n")

        # Create INCAR with NCORE = 4
        incar = os.path.join(temp_dir, "INCAR")
        with open(incar, "w") as f:
            f.write("NCORE = 4\n")

        nprocs = get_optimal_nprocs(temp_dir)

        # nprocs should be divisible by NCORE
        assert nprocs % 4 == 0 or nprocs < 4


class TestLocalRunner:
    """Test LocalRunner."""

    def test_init_defaults(self, monkeypatch):
        """Test default initialization."""
        # Clear VASP_COMMAND env var to test true default
        monkeypatch.delenv("VASP_COMMAND", raising=False)
        monkeypatch.delenv("VASP_EXECUTABLE", raising=False)
        runner = LocalRunner()

        assert runner.vasp_executable == "vasp_std"
        assert runner.nprocs == "auto"
        assert runner.mpi_command == "mpirun"
        assert not runner.background

    def test_init_from_env(self, monkeypatch):
        """Test initialization from VASP_EXECUTABLE environment variable."""
        monkeypatch.delenv("VASP_COMMAND", raising=False)
        monkeypatch.setenv("VASP_EXECUTABLE", "vasp_ncl")
        runner = LocalRunner()

        assert runner.vasp_executable == "vasp_ncl"

    def test_init_with_nprocs(self, monkeypatch):
        """Test initialization with explicit nprocs."""
        monkeypatch.delenv("VASP_COMMAND", raising=False)
        runner = LocalRunner(vasp_executable="vasp_gam", nprocs=4)

        assert runner.vasp_executable == "vasp_gam"
        assert runner.nprocs == 4

    def test_build_command(self, temp_dir, monkeypatch):
        """Test command building."""
        monkeypatch.delenv("VASP_COMMAND", raising=False)

        # Serial (nprocs=1, no MPI)
        runner = LocalRunner(vasp_executable="vasp_std", nprocs=1, mpi_command=None)
        assert runner._build_command(temp_dir) == "vasp_std"

        # With MPI
        runner = LocalRunner(vasp_executable="vasp_std", nprocs=4)
        assert runner._build_command(temp_dir) == "mpirun -np 4 vasp_std"

        # With extra args
        runner = LocalRunner(vasp_executable="vasp_std", nprocs=4, mpi_extra_args="--bind-to core")
        assert runner._build_command(temp_dir) == "mpirun -np 4 --bind-to core vasp_std"

    def test_verify_inputs_missing(self, temp_dir):
        """Test input verification with missing files."""
        runner = LocalRunner()

        with pytest.raises(VaspSetupError) as exc_info:
            runner._verify_inputs(temp_dir)

        assert "INCAR" in str(exc_info.value)

    def test_verify_inputs_complete(self, calc_dir):
        """Test input verification with complete files."""
        runner = LocalRunner()
        # Should not raise
        runner._verify_inputs(calc_dir)

    def test_status_not_started(self, temp_dir):
        """Test status for not started calculation."""
        runner = LocalRunner()

        status = runner.status(temp_dir)
        assert status.state == JobState.NOT_STARTED

    def test_status_complete(self, complete_calc_dir):
        """Test status for completed calculation."""
        runner = LocalRunner()

        status = runner.status(complete_calc_dir)
        assert status.state == JobState.COMPLETE

    def test_check_outcar_complete(self, complete_calc_dir):
        """Test OUTCAR completion check."""
        runner = LocalRunner()

        assert runner._check_outcar_complete(complete_calc_dir)

    def test_check_outcar_error(self, temp_dir):
        """Test OUTCAR error detection."""
        runner = LocalRunner()

        # Create OUTCAR with error
        with open(os.path.join(temp_dir, "OUTCAR"), "w") as f:
            f.write("VERY BAD NEWS! Something went wrong\n")

        error = runner._check_outcar_error(temp_dir)
        assert error == "VERY BAD NEWS!"

    def test_repr(self, monkeypatch):
        """Test string representation."""
        monkeypatch.delenv("VASP_COMMAND", raising=False)
        runner = LocalRunner(vasp_executable="vasp_ncl", nprocs=8)

        repr_str = repr(runner)
        assert "LocalRunner" in repr_str
        assert "vasp_ncl" in repr_str
        assert "nprocs=8" in repr_str


class TestSlurmRunner:
    """Test SlurmRunner."""

    def test_init_defaults(self):
        """Test default initialization."""
        runner = SlurmRunner()

        assert runner.nodes == 1
        assert runner.ntasks_per_node == 24  # Actual default
        assert runner.time == "24:00:00"
        assert "vasp_std" in runner.vasp_command

    def test_init_custom(self):
        """Test custom initialization."""
        runner = SlurmRunner(
            nodes=4,
            ntasks_per_node=32,
            partition="gpu",
            account="my-project",
            time="48:00:00",
            vasp_command="vasp_gpu",
        )

        assert runner.nodes == 4
        assert runner.ntasks_per_node == 32
        assert runner.partition == "gpu"
        assert runner.account == "my-project"

    def test_generate_script(self, temp_dir):
        """Test SLURM script generation."""
        runner = SlurmRunner(
            nodes=2,
            partition="normal",
        )

        # Method is named _create_script, not _generate_script
        script = runner._create_script(temp_dir)

        assert "#!/bin/bash" in script
        assert "#SBATCH --nodes=2" in script
        assert "#SBATCH --partition=normal" in script
        assert "#SBATCH --job-name=" in script  # job_name derived from directory
        assert "srun" in script or "vasp" in script

    def test_generate_script_with_modules(self, temp_dir):
        """Test script generation with module loading."""
        runner = SlurmRunner(modules=["vasp/6.3.0", "intel/2021"])

        script = runner._create_script(temp_dir)

        assert "module load vasp/6.3.0" in script
        assert "module load intel/2021" in script

    def test_status_not_started(self, temp_dir):
        """Test status for not started job."""
        runner = SlurmRunner()

        status = runner.status(temp_dir)
        assert status.state == JobState.NOT_STARTED

    def test_status_with_jobid_file(self, temp_dir, complete_calc_dir):
        """Test status with job ID file but no queue entry."""
        runner = SlurmRunner()

        # Create job ID file
        with open(os.path.join(complete_calc_dir, ".slurm_jobid"), "w") as f:
            f.write("12345")

        status = runner.status(complete_calc_dir)
        # Should check squeue, then fall back to OUTCAR
        assert status.state == JobState.COMPLETE

    def test_repr(self):
        """Test string representation."""
        runner = SlurmRunner(partition="gpu", nodes=4)

        repr_str = repr(runner)
        assert "SlurmRunner" in repr_str
        assert "gpu" in repr_str
        assert "4" in repr_str


class TestKubernetesRunner:
    """Test KubernetesRunner.

    Note: These tests require the kubernetes Python package.
    They will be skipped if the package is not installed.
    """

    @pytest.fixture(autouse=True)
    def check_kubernetes(self):
        """Skip tests if kubernetes package or config not available."""
        try:
            import kubernetes  # noqa: F401
            from kubernetes import config

            config.load_kube_config()
        except ImportError:
            pytest.skip("kubernetes package not installed")
        except kubernetes.config.config_exception.ConfigException:
            pytest.skip("kubernetes config not available")

    def test_init_defaults(self):
        """Test default initialization."""
        runner = KubernetesRunner(pvc_name="vasp-data")

        assert runner.pvc_name == "vasp-data"
        assert runner.namespace == "default"
        # attribute is 'image', not 'vasp_image'
        assert runner.image == "vasp:latest"

    def test_init_custom(self):
        """Test custom initialization."""
        runner = KubernetesRunner(
            pvc_name="my-pvc",
            namespace="compute",
            image="registry/vasp:6.3.0",  # 'image' not 'vasp_image'
            cpu_request="4",
            memory_request="16Gi",
            gpu_limit=1,  # 'gpu_limit' not 'gpu_request'
        )

        assert runner.pvc_name == "my-pvc"
        assert runner.namespace == "compute"
        assert runner.gpu_limit == 1

    def test_generate_job_name(self, temp_dir):
        """Test Kubernetes Job name generation."""
        runner = KubernetesRunner(pvc_name="vasp-data")

        job_name = runner._generate_job_name(temp_dir)

        # Job name should start with 'vasp-' and be K8s compatible
        assert job_name.startswith("vasp-")
        assert all(c.isalnum() or c == "-" for c in job_name)

    def test_status_not_started(self, temp_dir):
        """Test status for not started job."""
        runner = KubernetesRunner(pvc_name="vasp-data")

        status = runner.status(temp_dir)
        assert status.state == JobState.NOT_STARTED

    def test_repr(self):
        """Test string representation."""
        runner = KubernetesRunner(pvc_name="my-pvc", namespace="compute")

        repr_str = repr(runner)
        assert "KubernetesRunner" in repr_str
        # Check for namespace in repr
        assert "compute" in repr_str


class TestRunnerBase:
    """Test base Runner functionality."""

    def test_check_outcar_complete_missing(self, temp_dir):
        """Test completion check with missing OUTCAR."""
        runner = MockRunner(results=MockResults(energy=-10.0))

        assert not runner._check_outcar_complete(temp_dir)

    def test_check_outcar_complete_incomplete(self, temp_dir):
        """Test completion check with incomplete OUTCAR."""
        runner = MockRunner(results=MockResults(energy=-10.0))

        # Write incomplete OUTCAR
        with open(os.path.join(temp_dir, "OUTCAR"), "w") as f:
            f.write("Some output but not complete\n")

        assert not runner._check_outcar_complete(temp_dir)

    def test_get_logs(self, complete_calc_dir):
        """Test getting logs from OUTCAR."""
        runner = MockRunner(results=MockResults(energy=-10.0))

        logs = runner.get_logs(complete_calc_dir, tail_lines=10)

        assert isinstance(logs, str)
        assert len(logs) > 0

    def test_get_logs_missing(self, temp_dir):
        """Test getting logs with no OUTCAR."""
        runner = MockRunner(results=MockResults(energy=-10.0))

        logs = runner.get_logs(temp_dir)

        assert "No OUTCAR found" in logs
