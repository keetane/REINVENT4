from __future__ import annotations

import os
import signal
import socket
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.models.job import JobStatus
from app.reinvent_adapter import AdapterConfigurationError, build_reinvent_command
from app.services.naming import confined_path, result_filename
from app.services.results import ResultValidationError, validate_and_finalize_csv
from app.services.workspaces import relative_to_repository


def worker_identity() -> str:
    return f"{socket.gethostname()}-{os.getpid()}-{str(uuid4())[:8]}"


def _allowed_environment() -> dict[str, str]:
    allowed = {"PATH", "LANG", "LC_ALL", "TMPDIR", "SYSTEMROOT"}
    environment = {key: value for key, value in os.environ.items() if key in allowed}
    environment["PYTHONUNBUFFERED"] = "1"
    return environment


@dataclass(frozen=True, slots=True)
class FakeCommand:
    argv: list[str]
    raw_output: Path
    process_log: Path
    public_log: Path
    environment: dict[str, str]
    config_path: Path | None = None
    model_path: str | None = None


def fake_command(job: JobRecord, job_dir: Path) -> FakeCommand:
    conditions = job.run_conditions
    delay = max(0.0, min(float(conditions.get("fake_delay_seconds", 0.05)), 60.0))
    rows = max(1, min(int(conditions.get("fake_rows", 2)), 1_000))
    exit_code = max(0, min(int(conditions.get("fake_exit_code", 0)), 255))
    output_mode = str(conditions.get("fake_output_mode", "valid"))
    if output_mode not in {"valid", "missing", "empty", "missing-smiles"}:
        raise ValueError("fake_output_mode is not supported")
    raw_output = job_dir / "output" / "raw.csv"
    script = Path(__file__).with_name("fake_job.py").resolve()
    log_path = job_dir / "logs" / "worker.log"
    return FakeCommand(
        argv=[
            sys.executable,
            str(script),
            "--output",
            str(raw_output),
            "--delay",
            str(delay),
            "--rows",
            str(rows),
            "--exit-code",
            str(exit_code),
            "--output-mode",
            output_mode,
        ],
        raw_output=raw_output,
        process_log=log_path,
        public_log=log_path,
        environment={},
    )


def execution_command(settings: Settings, job: JobRecord, job_dir: Path):
    if settings.worker_runner == "fake":
        return fake_command(job, job_dir)
    if settings.worker_runner == "reinvent":
        return build_reinvent_command(settings, job, job_dir)
    raise ValueError("configured worker runner is not supported")


class Worker:
    def __init__(self, settings: Settings, *, worker_id: str | None = None):
        self.settings = settings
        self.repository = JobRepository(settings)
        self.worker_id = worker_id or worker_identity()

    def run_once(self) -> bool:
        job = self.repository.claim_next(self.worker_id)
        if job is None:
            return False
        self._execute(job)
        return True

    def run_forever(self) -> None:
        while True:
            if not self.run_once():
                reconcile_active_jobs(self.settings)
                time.sleep(self.settings.worker_poll_interval_seconds)

    def _execute(self, job: JobRecord) -> None:
        exit_code: int | None = None
        try:
            job_dir = confined_path(self.settings.repository_root, job.job_dir)
            command = execution_command(self.settings, job, job_dir)
            environment = _allowed_environment()
            environment.update(command.environment)
            with command.process_log.open("ab") as log_handle:
                process = subprocess.Popen(
                    command.argv,
                    cwd=job_dir,
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=log_handle,
                    stderr=subprocess.STDOUT,
                    shell=False,
                    start_new_session=True,
                )
                self.repository.set_process(
                    job.id,
                    process_id=process.pid,
                    log_path=relative_to_repository(self.settings, command.public_log),
                    command=command.argv,
                    config_path=(
                        relative_to_repository(self.settings, command.config_path)
                        if command.config_path
                        else None
                    ),
                    model_path=command.model_path,
                )
                exit_code = self._wait(job.id, process)
            current = self.repository.get(job.id)
            if current.status == JobStatus.CANCELLING:
                self.repository.finish(
                    job.id, target=JobStatus.CANCELLED, exit_code=exit_code
                )
                return
            if exit_code != 0:
                self.repository.finish(
                    job.id,
                    target=JobStatus.FAILED,
                    exit_code=exit_code,
                    error_summary=f"The worker command exited with code {exit_code}.",
                    error_kind="nonzero_exit",
                )
                return
            self.repository.set_progress(job.id, "validating output")
            created_at = datetime.fromisoformat(job.created_at)
            filename = result_filename(
                comment_slug=job.comment_slug,
                created_at=created_at,
                timezone_name=self.settings.app_timezone,
            )
            final_path = job_dir / "output" / filename
            validate_and_finalize_csv(command.raw_output, final_path, job_dir)
            self.repository.finish(
                job.id,
                target=JobStatus.COMPLETED,
                exit_code=exit_code,
                result_csv_path=relative_to_repository(self.settings, final_path),
            )
        except ResultValidationError as error:
            self._fail_if_active(job.id, exit_code, str(error), "invalid_output")
        except AdapterConfigurationError as error:
            self._fail_if_active(
                job.id, exit_code, str(error), "adapter_configuration"
            )
        except PermissionError:
            self._fail_if_active(
                job.id,
                exit_code,
                "The worker does not have permission to access a required path.",
                "permission_denied",
            )
        except FileNotFoundError:
            self._fail_if_active(
                job.id,
                exit_code,
                "A required executable or file is missing.",
                "required_file_missing",
            )
        except Exception as error:
            self._fail_if_active(
                job.id,
                exit_code,
                "The worker could not complete this job.",
                type(error).__name__,
            )

    def _wait(self, job_id, process: subprocess.Popen) -> int:
        while True:
            exit_code = process.poll()
            if exit_code is not None:
                return exit_code
            if self.repository.get(job_id).status == JobStatus.CANCELLING:
                self._terminate_process_group(process)
                return process.wait()
            time.sleep(0.05)

    def _terminate_process_group(self, process: subprocess.Popen) -> None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            return
        try:
            process.wait(timeout=self.settings.cancel_grace_seconds)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

    def _fail_if_active(
        self,
        job_id,
        exit_code: int | None,
        summary: str,
        kind: str,
    ) -> None:
        current = self.repository.get(job_id)
        if current.status == JobStatus.CANCELLING:
            self.repository.finish(
                job_id, target=JobStatus.CANCELLED, exit_code=exit_code
            )
        elif current.status == JobStatus.RUNNING:
            self.repository.finish(
                job_id,
                target=JobStatus.FAILED,
                exit_code=exit_code,
                error_summary=summary,
                error_kind=kind,
            )


def process_is_alive(process_id: int | None) -> bool:
    if process_id is None:
        return False
    try:
        os.kill(process_id, 0)
    except (OSError, ValueError):
        return False
    return True


def reconcile_active_jobs(settings: Settings) -> list[str]:
    repository = JobRepository(settings)
    interrupted: list[str] = []
    now = datetime.now(timezone.utc)
    for job in repository.active_jobs():
        if process_is_alive(job.process_id):
            continue
        updated_at = datetime.fromisoformat(job.updated_at)
        if (now - updated_at).total_seconds() < settings.reconciliation_grace_seconds:
            continue
        repository.interrupt(
            job.id,
            "The previous worker is no longer running; the job was not restarted automatically.",
        )
        interrupted.append(str(job.id))
    return interrupted
