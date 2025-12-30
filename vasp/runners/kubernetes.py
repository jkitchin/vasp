"""Kubernetes runner for VASP execution on K8s clusters.

Submits VASP jobs as Kubernetes Job resources. Requires:
- kubernetes Python package: pip install kubernetes
- Configured kubectl context or in-cluster config
- PersistentVolumeClaim for calculation directories
"""

from __future__ import annotations

import os
from datetime import datetime
from typing import TYPE_CHECKING, Any

from ..exceptions import VaspQueued, VaspRunning, VaspSubmitted
from .base import JobState, JobStatus, Runner

if TYPE_CHECKING:
    pass

# Check for kubernetes package
try:
    from kubernetes import client, config
    from kubernetes.client.rest import ApiException

    HAS_KUBERNETES = True
except ImportError:
    HAS_KUBERNETES = False
    client = None
    config = None
    ApiException = Exception


class KubernetesRunner(Runner):
    """Run VASP on Kubernetes cluster.

    Submits VASP calculations as Kubernetes Job resources.
    Input/output files must be on shared storage (PVC).

    Args:
        namespace: Kubernetes namespace.
        image: Container image with VASP.
        pvc_name: PersistentVolumeClaim name for data.
        pvc_mount_path: Mount path in container.
        cpu_request: CPU request (e.g., '4').
        cpu_limit: CPU limit (defaults to request).
        memory_request: Memory request (e.g., '16Gi').
        memory_limit: Memory limit (defaults to request).
        gpu_limit: Number of GPUs (for vasp_gpu).
        node_selector: Dict of node labels for scheduling.
        tolerations: List of tolerations for scheduling.
        service_account: ServiceAccount name.
        image_pull_secrets: List of image pull secret names.
        vasp_command: Command to run VASP.
        env_vars: Environment variables dict.
        ttl_seconds_after_finished: Job cleanup time (default 1 day).
        active_deadline_seconds: Job timeout.
        backoff_limit: Retry count (default 0 = no retries).
        use_in_cluster_config: Use in-cluster config (for pods).

    Example:
        >>> runner = KubernetesRunner(
        ...     namespace='vasp-jobs',
        ...     image='my-registry/vasp:6.3.0',
        ...     pvc_name='vasp-data',
        ...     cpu_request='8',
        ...     memory_request='32Gi',
        ... )
        >>>
        >>> calc = Vasp('/data/calculations/my_calc', runner=runner, ...)
    """

    def __init__(
        self,
        namespace: str = "default",
        image: str = "vasp:latest",
        pvc_name: str = "vasp-pvc",
        pvc_mount_path: str = "/data",
        cpu_request: str = "4",
        cpu_limit: str | None = None,
        memory_request: str = "8Gi",
        memory_limit: str | None = None,
        gpu_limit: int | None = None,
        node_selector: dict[str, str] | None = None,
        tolerations: list[dict] | None = None,
        service_account: str | None = None,
        image_pull_secrets: list[str] | None = None,
        vasp_command: str = "mpirun -np $NPROCS vasp_std",
        env_vars: dict[str, str] | None = None,
        ttl_seconds_after_finished: int = 86400,
        active_deadline_seconds: int | None = None,
        backoff_limit: int = 0,
        use_in_cluster_config: bool = False,
    ):
        if not HAS_KUBERNETES:
            raise ImportError("kubernetes package required: pip install kubernetes")

        self.namespace = namespace
        self.image = image
        self.pvc_name = pvc_name
        self.pvc_mount_path = pvc_mount_path
        self.cpu_request = cpu_request
        self.cpu_limit = cpu_limit or cpu_request
        self.memory_request = memory_request
        self.memory_limit = memory_limit or memory_request
        self.gpu_limit = gpu_limit
        self.node_selector = node_selector or {}
        self.tolerations = tolerations or []
        self.service_account = service_account
        self.image_pull_secrets = image_pull_secrets or []
        self.vasp_command = vasp_command
        self.env_vars = env_vars or {}
        self.ttl_seconds_after_finished = ttl_seconds_after_finished
        self.active_deadline_seconds = active_deadline_seconds
        self.backoff_limit = backoff_limit

        # Initialize Kubernetes client
        if use_in_cluster_config:
            config.load_incluster_config()
        else:
            config.load_kube_config()

        self.batch_v1 = client.BatchV1Api()
        self.core_v1 = client.CoreV1Api()

    def run(self, directory: str) -> JobStatus:
        """Submit VASP job to Kubernetes."""
        current = self.status(directory)

        if current.state == JobState.COMPLETE:
            return current
        if current.state == JobState.QUEUED:
            raise VaspQueued(message=f"Pod pending: {current.jobid}", jobid=current.jobid)
        if current.state == JobState.RUNNING:
            raise VaspRunning(message=f"Pod running: {current.jobid}", jobid=current.jobid)
        if current.state == JobState.FAILED:
            self._cleanup_old_job(directory)

        # Submit new job
        job_name = self._submit(directory)
        raise VaspSubmitted(jobid=job_name)

    def status(self, directory: str) -> JobStatus:
        """Check Kubernetes job status."""
        job_name = self._read_job_name(directory)

        if job_name:
            try:
                job = self.batch_v1.read_namespaced_job(name=job_name, namespace=self.namespace)
                return self._parse_job_status(job, directory)
            except ApiException as e:
                if e.status == 404:
                    pass  # Job deleted, check output files
                else:
                    raise

        return self._check_output_files(directory)

    def cancel(self, directory: str) -> bool:
        """Delete Kubernetes job."""
        job_name = self._read_job_name(directory)
        if not job_name:
            return True

        try:
            self.batch_v1.delete_namespaced_job(
                name=job_name,
                namespace=self.namespace,
                body=client.V1DeleteOptions(propagation_policy="Foreground"),
            )
            return True
        except ApiException:
            return False

    def get_logs(self, directory: str, tail_lines: int = 100) -> str:
        """Get logs from VASP pod."""
        job_name = self._read_job_name(directory)
        if not job_name:
            return super().get_logs(directory, tail_lines)

        try:
            pods = self.core_v1.list_namespaced_pod(
                namespace=self.namespace, label_selector=f"job-name={job_name}"
            )

            if not pods.items:
                return "No pods found for job"

            pod = pods.items[0]
            return str(
                self.core_v1.read_namespaced_pod_log(
                    name=pod.metadata.name, namespace=self.namespace, tail_lines=tail_lines
                )
            )
        except ApiException as e:
            return f"Error getting logs: {e}"

    def _submit(self, directory: str) -> str:
        """Create and submit Kubernetes Job."""
        job_name = self._generate_job_name(directory)
        work_dir = self._get_work_dir(directory)

        job = self._create_job_manifest(job_name, work_dir)

        self.batch_v1.create_namespaced_job(namespace=self.namespace, body=job)

        self._write_job_name(directory, job_name)
        return job_name

    def _create_job_manifest(self, job_name: str, work_dir: str) -> Any:
        """Create Kubernetes Job manifest."""
        env = [client.V1EnvVar(name=k, value=v) for k, v in self.env_vars.items()]
        env.append(client.V1EnvVar(name="NPROCS", value=self.cpu_request))

        resources = client.V1ResourceRequirements(
            requests={"cpu": self.cpu_request, "memory": self.memory_request},
            limits={"cpu": self.cpu_limit, "memory": self.memory_limit},
        )

        if self.gpu_limit:
            resources.limits["nvidia.com/gpu"] = str(self.gpu_limit)

        volume_mount = client.V1VolumeMount(name="vasp-data", mount_path=self.pvc_mount_path)

        container = client.V1Container(
            name="vasp",
            image=self.image,
            command=["/bin/bash", "-c"],
            args=[f"cd {self.pvc_mount_path}/{work_dir} && {self.vasp_command}"],
            env=env,
            resources=resources,
            volume_mounts=[volume_mount],
            working_dir=f"{self.pvc_mount_path}/{work_dir}",
        )

        volume = client.V1Volume(
            name="vasp-data",
            persistent_volume_claim=client.V1PersistentVolumeClaimVolumeSource(
                claim_name=self.pvc_name
            ),
        )

        pull_secrets = [
            client.V1LocalObjectReference(name=s) for s in self.image_pull_secrets
        ] or None

        pod_spec = client.V1PodSpec(
            containers=[container],
            volumes=[volume],
            restart_policy="Never",
            node_selector=self.node_selector or None,
            tolerations=[client.V1Toleration(**t) for t in self.tolerations] or None,
            service_account_name=self.service_account,
            image_pull_secrets=pull_secrets,
        )

        job_spec = client.V1JobSpec(
            template=client.V1PodTemplateSpec(
                metadata=client.V1ObjectMeta(labels={"app": "vasp", "job-name": job_name}),
                spec=pod_spec,
            ),
            backoff_limit=self.backoff_limit,
            ttl_seconds_after_finished=self.ttl_seconds_after_finished,
            active_deadline_seconds=self.active_deadline_seconds,
        )

        return client.V1Job(
            api_version="batch/v1",
            kind="Job",
            metadata=client.V1ObjectMeta(
                name=job_name, namespace=self.namespace, labels={"app": "vasp"}
            ),
            spec=job_spec,
        )

    def _parse_job_status(self, job: Any, directory: str) -> JobStatus:
        """Parse Kubernetes Job status."""
        status = job.status
        job_name = job.metadata.name

        if status.succeeded and status.succeeded > 0:
            return JobStatus(JobState.COMPLETE, jobid=job_name)

        if status.failed and status.failed > 0:
            return JobStatus(JobState.FAILED, jobid=job_name, message=self._get_failure_reason(job))

        if status.active and status.active > 0:
            pods = self.core_v1.list_namespaced_pod(
                namespace=self.namespace, label_selector=f"job-name={job_name}"
            )

            if pods.items:
                pod_phase = pods.items[0].status.phase
                if pod_phase == "Pending":
                    return JobStatus(JobState.QUEUED, jobid=job_name)
                elif pod_phase == "Running":
                    return JobStatus(JobState.RUNNING, jobid=job_name)

        return JobStatus(JobState.QUEUED, jobid=job_name)

    def _get_failure_reason(self, job: Any) -> str:
        """Extract failure reason from job."""
        if job.status.conditions:
            for condition in job.status.conditions:
                if condition.type == "Failed":
                    return condition.message or condition.reason or "Unknown"
        return "Job failed"

    def _check_output_files(self, directory: str) -> JobStatus:
        """Check calculation status from output files."""
        if self._check_outcar_complete(directory):
            return JobStatus(JobState.COMPLETE)

        error = self._check_outcar_error(directory)
        if error:
            return JobStatus(JobState.FAILED, message=error)

        outcar = os.path.join(directory, "OUTCAR")
        if os.path.exists(outcar):
            return JobStatus(JobState.FAILED, message="OUTCAR incomplete")

        return JobStatus(JobState.NOT_STARTED)

    def _generate_job_name(self, directory: str) -> str:
        """Generate unique K8s-compatible job name."""
        base = os.path.basename(os.path.abspath(directory)).lower()
        base = "".join(c if c.isalnum() or c == "-" else "-" for c in base)
        base = base.strip("-")[:40]
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        return f"vasp-{base}-{timestamp}"

    def _get_work_dir(self, directory: str) -> str:
        """Get directory relative to PVC mount."""
        return os.path.abspath(directory).lstrip("/")

    def _read_job_name(self, directory: str) -> str | None:
        """Read K8s job name from tracking file."""
        path = os.path.join(directory, ".k8s_job_name")
        if os.path.exists(path):
            with open(path) as f:
                return f.read().strip()
        return None

    def _write_job_name(self, directory: str, job_name: str) -> None:
        """Save K8s job name to tracking file."""
        path = os.path.join(directory, ".k8s_job_name")
        with open(path, "w") as f:
            f.write(job_name)

    def _cleanup_old_job(self, directory: str) -> None:
        """Delete old job before resubmission."""
        job_name = self._read_job_name(directory)
        if job_name:
            try:
                self.batch_v1.delete_namespaced_job(
                    name=job_name,
                    namespace=self.namespace,
                    body=client.V1DeleteOptions(propagation_policy="Background"),
                )
            except ApiException:
                pass

            path = os.path.join(directory, ".k8s_job_name")
            if os.path.exists(path):
                os.remove(path)

    def __repr__(self) -> str:
        return f"KubernetesRunner(namespace={self.namespace!r}, " f"image={self.image!r})"
