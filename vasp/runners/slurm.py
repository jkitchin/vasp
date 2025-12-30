"""SLURM runner for VASP execution on HPC clusters.

Submits VASP jobs to SLURM scheduler and monitors their status.
Designed for non-blocking async operation.
"""

from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

from ..exceptions import VaspQueued, VaspRunning, VaspSubmitted
from .base import JobState, JobStatus, Runner

if TYPE_CHECKING:
    pass


class SlurmRunner(Runner):
    """Run VASP via SLURM scheduler.

    Submits jobs to SLURM and monitors status. All operations are
    non-blocking - run() submits and returns immediately.

    Args:
        partition: SLURM partition name.
        nodes: Number of nodes.
        ntasks_per_node: MPI tasks per node.
        time: Wall time limit (HH:MM:SS format).
        memory: Memory per node (e.g., '64G').
        account: SLURM account for billing.
        qos: Quality of service.
        vasp_command: Command to run VASP. Defaults to 'srun $VASP_COMMAND'
            if VASP_COMMAND is set, or 'srun vasp_std' otherwise.
        modules: List of modules to load before running.
        extra_sbatch: Additional #SBATCH directives as list of strings.
        constraint: Node constraint (e.g., 'gpu' or 'skylake').

    Example:
        >>> runner = SlurmRunner(
        ...     partition='compute',
        ...     nodes=2,
        ...     ntasks_per_node=48,
        ...     time='24:00:00',
        ...     modules=['vasp/6.3.0'],
        ... )
        >>>
        >>> calc = Vasp('my_calc', runner=runner, atoms=atoms)
        >>> try:
        ...     energy = calc.potential_energy
        ... except VaspSubmitted as e:
        ...     print(f"Submitted job {e.jobid}")
    """

    def __init__(
        self,
        partition: str = 'normal',
        nodes: int = 1,
        ntasks_per_node: int = 24,
        time: str = '24:00:00',
        memory: str | None = None,
        account: str | None = None,
        qos: str | None = None,
        vasp_command: str | None = None,
        modules: list[str] | None = None,
        extra_sbatch: list[str] | None = None,
        constraint: str | None = None,
    ):
        self.partition = partition
        self.nodes = nodes
        self.ntasks_per_node = ntasks_per_node
        self.time = time
        self.memory = memory
        self.account = account
        self.qos = qos
        if vasp_command is not None:
            self.vasp_command = vasp_command
        else:
            base_cmd = os.environ.get('VASP_COMMAND', 'vasp_std')
            self.vasp_command = f'srun {base_cmd}'
        self.modules = modules or []
        self.extra_sbatch = extra_sbatch or []
        self.constraint = constraint

    def run(self, directory: str) -> JobStatus:
        """Submit VASP job to SLURM."""
        current = self.status(directory)

        if current.state == JobState.COMPLETE:
            return current
        if current.state == JobState.QUEUED:
            raise VaspQueued(message=f"Job {current.jobid} queued", jobid=current.jobid)
        if current.state == JobState.RUNNING:
            raise VaspRunning(message=f"Job {current.jobid} running", jobid=current.jobid)
        if current.state == JobState.FAILED:
            # Allow resubmission - clean up old job file
            self._cleanup_old_job(directory)

        # Submit new job
        jobid = self._submit(directory)
        raise VaspSubmitted(jobid=jobid)

    def status(self, directory: str) -> JobStatus:
        """Check SLURM job status."""
        jobid = self._read_jobid(directory)

        if jobid:
            state = self._query_slurm(jobid)
            if state:
                return JobStatus(state, jobid=jobid)

        # Job not in SLURM - check output files
        return self._check_output_files(directory)

    def cancel(self, directory: str) -> bool:
        """Cancel SLURM job."""
        jobid = self._read_jobid(directory)
        if not jobid:
            return True

        try:
            result = subprocess.run(
                ['scancel', jobid],
                capture_output=True,
                text=True,
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def get_logs(self, directory: str, tail_lines: int = 100) -> str:
        """Get SLURM job output."""
        # Try SLURM output file first
        jobid = self._read_jobid(directory)
        if jobid:
            slurm_out = os.path.join(directory, f'slurm-{jobid}.out')
            if os.path.exists(slurm_out):
                with open(slurm_out) as f:
                    lines = f.readlines()
                    return ''.join(lines[-tail_lines:])

        # Fall back to OUTCAR
        return super().get_logs(directory, tail_lines)

    def _submit(self, directory: str) -> str:
        """Create and submit SLURM job script."""
        script = self._create_script(directory)
        script_path = os.path.join(directory, 'submit.slurm')

        with open(script_path, 'w') as f:
            f.write(script)

        result = subprocess.run(
            ['sbatch', script_path],
            cwd=directory,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(f"sbatch failed: {result.stderr}")

        # Parse "Submitted batch job 123456"
        jobid = result.stdout.strip().split()[-1]
        self._write_jobid(directory, jobid)

        return jobid

    def _query_slurm(self, jobid: str) -> JobState | None:
        """Query SLURM for job state."""
        try:
            result = subprocess.run(
                ['squeue', '-j', jobid, '-h', '-o', '%T'],
                capture_output=True,
                text=True,
                timeout=30,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None

        if result.returncode != 0:
            return None

        state = result.stdout.strip()
        mapping = {
            'PENDING': JobState.QUEUED,
            'CONFIGURING': JobState.QUEUED,
            'RUNNING': JobState.RUNNING,
            'COMPLETING': JobState.RUNNING,
            'SUSPENDED': JobState.QUEUED,
        }
        return mapping.get(state)

    def _check_output_files(self, directory: str) -> JobStatus:
        """Check calculation status from output files."""
        if self._check_outcar_complete(directory):
            return JobStatus(JobState.COMPLETE)

        error = self._check_outcar_error(directory)
        if error:
            return JobStatus(JobState.FAILED, message=error)

        outcar = os.path.join(directory, 'OUTCAR')
        if os.path.exists(outcar):
            return JobStatus(JobState.FAILED, message="OUTCAR incomplete")

        return JobStatus(JobState.NOT_STARTED)

    def _create_script(self, directory: str) -> str:
        """Generate SLURM batch script."""
        job_name = os.path.basename(os.path.abspath(directory))[:50]

        lines = [
            '#!/bin/bash',
            f'#SBATCH --job-name={job_name}',
            f'#SBATCH --partition={self.partition}',
            f'#SBATCH --nodes={self.nodes}',
            f'#SBATCH --ntasks-per-node={self.ntasks_per_node}',
            f'#SBATCH --time={self.time}',
            '#SBATCH --output=slurm-%j.out',
            '#SBATCH --error=slurm-%j.err',
        ]

        if self.memory:
            lines.append(f'#SBATCH --mem={self.memory}')
        if self.account:
            lines.append(f'#SBATCH --account={self.account}')
        if self.qos:
            lines.append(f'#SBATCH --qos={self.qos}')
        if self.constraint:
            lines.append(f'#SBATCH --constraint={self.constraint}')

        for directive in self.extra_sbatch:
            lines.append(f'#SBATCH {directive}')

        lines.append('')
        lines.append('# Load modules')
        for mod in self.modules:
            lines.append(f'module load {mod}')

        lines.append('')
        lines.append('# Run VASP')
        lines.append(self.vasp_command)

        return '\n'.join(lines) + '\n'

    def _read_jobid(self, directory: str) -> str | None:
        """Read SLURM job ID from tracking file."""
        path = os.path.join(directory, '.slurm_jobid')
        if os.path.exists(path):
            with open(path) as f:
                return f.read().strip()
        return None

    def _write_jobid(self, directory: str, jobid: str) -> None:
        """Save SLURM job ID to tracking file."""
        path = os.path.join(directory, '.slurm_jobid')
        with open(path, 'w') as f:
            f.write(jobid)

    def _cleanup_old_job(self, directory: str) -> None:
        """Remove old job tracking file before resubmission."""
        path = os.path.join(directory, '.slurm_jobid')
        if os.path.exists(path):
            os.remove(path)

    def __repr__(self) -> str:
        return (
            f"SlurmRunner(partition={self.partition!r}, "
            f"nodes={self.nodes}, time={self.time!r})"
        )
