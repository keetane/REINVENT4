from __future__ import annotations

import json
import csv
import os
import re
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
from app.db.models import ModelRepository
from app.db.qsar_models import QsarModelRepository
from app.models.job import JobStatus
from app.reinvent_adapter import (
    AdapterConfigurationError,
    build_reinvent_command,
    build_rl_command,
    build_transfer_learning_command,
)
from app.services.naming import confined_path, result_filename
from app.services.results import ResultValidationError, validate_and_finalize_csv
from app.services.workspaces import relative_to_repository
from app.services.smiles_standardization import standardize_smiles


def _write_rl_audit(summary_path: Path) -> Path:
    audit_path = summary_path.with_name("rl-standardization-audit.csv")
    with summary_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        source_rows = list(reader)
        source_columns = list(reader.fieldnames or [])
    extra = [
        "Generated_SMILES", "Standardized_SMILES", "MolecularWeight", "FormalCharge",
        "QSAR_Prediction", "AD_Status", "QSAR_Score", "Total_Score", "Exclusion_Reason",
    ]
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=extra + source_columns)
        writer.writeheader()
        for source in source_rows:
            generated = source.get("SMILES") or source.get("smiles") or ""
            result = standardize_smiles(generated)
            raw_key = next((key for key in source if key.startswith("QSAR") and "raw" in key.lower()), None)
            score_key = next((key for key in source if key.startswith("QSAR:") and "raw" not in key.lower()), None)
            ad_key = next((key for key in source if key.startswith("QSAR AD:") and "raw" not in key.lower()), None)
            ad_value = source.get(ad_key, "") if ad_key else ""
            writer.writerow({
                "Generated_SMILES": generated,
                "Standardized_SMILES": result.standardized_smiles or "",
                "MolecularWeight": result.molecular_weight if result.molecular_weight is not None else "",
                "FormalCharge": result.formal_charge if result.formal_charge is not None else "",
                "QSAR_Prediction": source.get(raw_key, "") if raw_key else "",
                "AD_Status": "Excluded" if not result.accepted else "In" if str(ad_value) == "1.0" else "Out" if str(ad_value) == "0.0" else "Not evaluated",
                "QSAR_Score": source.get(score_key, "") if score_key else "",
                "Total_Score": source.get("Score", source.get("total_score", "")),
                "Exclusion_Reason": result.reason or "",
                **source,
            })
    return audit_path


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
    raw_output: Path | None
    process_log: Path
    public_log: Path
    environment: dict[str, str]
    config_path: Path | None = None
    model_path: str | None = None
    output_model_path: Path | None = None


def qsar_training_command(settings: Settings, job: JobRecord, job_dir: Path) -> FakeCommand:
    conditions = job.run_conditions
    output_model = confined_path(settings.repository_root, conditions["output_model_path"])
    config_path = confined_path(settings.repository_root, conditions["config_path"])
    log_path = job_dir / "logs" / "qsar-training.log"
    backend_root = Path(__file__).resolve().parents[1]
    module = "worker.chemprop_v2_training" if conditions["algorithm"] == "chemprop_v2" else "worker.qsartuna_training" if conditions["algorithm"] == "qsartuna" else "worker.qsar_training"
    return FakeCommand(
        argv=[
            str(settings.reinvent_runtime_python),
            "-m",
            module,
            "--config",
            str(config_path),
            "--repository-root",
            str(settings.repository_root),
        ],
        raw_output=None,
        process_log=log_path,
        public_log=log_path,
        environment={"PYTHONPATH": str(backend_root)},
        config_path=config_path,
        model_path=conditions["output_model_path"],
        output_model_path=output_model,
    )


def qsar_comparison_command(settings: Settings, job: JobRecord, job_dir: Path) -> FakeCommand:
    conditions = job.run_conditions
    config_path = confined_path(settings.repository_root, conditions["config_path"])
    log_path = job_dir / "logs" / "qsar-comparison.log"
    backend_root = Path(__file__).resolve().parents[1]
    return FakeCommand(
        argv=[
            str(settings.reinvent_runtime_python), "-m", "worker.qsar_comparison",
            "--config", str(config_path), "--repository-root", str(settings.repository_root),
        ],
        raw_output=None,
        process_log=log_path,
        public_log=log_path,
        environment={"PYTHONPATH": str(backend_root)},
        config_path=config_path,
    )


def fake_qsar_training_command(settings: Settings, job: JobRecord, job_dir: Path) -> FakeCommand:
    conditions = job.run_conditions
    output_model = confined_path(settings.repository_root, conditions["output_model_path"])
    evaluation = confined_path(settings.repository_root, conditions["evaluation_path"])
    predictions = confined_path(settings.repository_root, conditions["predictions_path"])
    config_path = confined_path(settings.repository_root, conditions["config_path"])
    log_path = job_dir / "logs" / "qsar-training.log"
    script = Path(__file__).with_name("fake_qsar_training.py").resolve()
    return FakeCommand(
        argv=[
            sys.executable,
            str(script),
            "--model",
            str(output_model),
            "--evaluation",
            str(evaluation),
            "--predictions",
            str(predictions),
            "--training-count",
            str(int(conditions["valid_count"])),
            "--config",
            str(config_path),
        ],
        raw_output=None,
        process_log=log_path,
        public_log=log_path,
        environment={},
        config_path=config_path,
        model_path=conditions["output_model_path"],
        output_model_path=output_model,
    )


def fake_training_command(settings: Settings, job: JobRecord, job_dir: Path) -> FakeCommand:
    output_model = (settings.repository_root / job.run_conditions["output_model_path"]).resolve()
    script = Path(__file__).with_name("fake_training.py").resolve()
    log_path = job_dir / "logs" / "transfer-learning.log"
    config_path = job_dir / "config" / "transfer-learning.toml"
    config_path.write_text(
        'run_type = "transfer_learning"\n[parameters]\n'
        f'num_epochs = {int(job.run_conditions.get("num_epochs", 3))}\n'
        f'output_model_file = "{output_model}"\n',
        encoding="utf-8",
    )
    return FakeCommand(
        argv=[
            sys.executable,
            str(script),
            "--output",
            str(output_model),
            "--epochs",
            str(int(job.run_conditions.get("num_epochs", 3))),
        ],
        raw_output=None,
        process_log=log_path,
        public_log=log_path,
        environment={},
        config_path=config_path,
        model_path=job.run_conditions["output_model_path"],
        output_model_path=output_model,
    )


def fake_rl_command(settings: Settings, job: JobRecord, job_dir: Path) -> FakeCommand:
    conditions = job.run_conditions
    output_model = confined_path(settings.repository_root, conditions["output_model_path"])
    summary = confined_path(settings.repository_root, conditions["summary_csv_path"])
    log_path = job_dir / "logs" / "reinforcement-learning.log"
    config_path = job_dir / "config" / "reinforcement-learning.toml"
    config_path.write_text('run_type = "staged_learning"\n', encoding="utf-8")
    script = Path(__file__).with_name("fake_rl.py").resolve()
    return FakeCommand(
        argv=[sys.executable, str(script), "--checkpoint", str(output_model), "--summary", str(summary)],
        raw_output=None,
        process_log=log_path,
        public_log=log_path,
        environment={},
        config_path=config_path,
        model_path=conditions["output_model_path"],
        output_model_path=output_model,
    )
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
        if job.job_type == "qsar_comparison":
            return qsar_comparison_command(settings, job, job_dir)
        if job.job_type == "reinforcement_learning":
            return fake_rl_command(settings, job, job_dir)
        if job.job_type == "qsar_training":
            return fake_qsar_training_command(settings, job, job_dir)
        if job.job_type == "transfer_learning":
            return fake_training_command(settings, job, job_dir)
        return fake_command(job, job_dir)
    if settings.worker_runner == "reinvent":
        if job.job_type == "qsar_comparison":
            return qsar_comparison_command(settings, job, job_dir)
        if job.job_type == "reinforcement_learning":
            return build_rl_command(settings, job, job_dir)
        if job.job_type == "qsar_training":
            return qsar_training_command(settings, job, job_dir)
        if job.job_type == "transfer_learning":
            return build_transfer_learning_command(settings, job, job_dir)
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
            visible_log = (
                command.process_log
                if job.job_type in {"transfer_learning", "qsar_training", "reinforcement_learning"}
                else command.public_log
            )
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
                    log_path=relative_to_repository(self.settings, visible_log),
                    command=command.argv,
                    config_path=(
                        relative_to_repository(self.settings, command.config_path)
                        if command.config_path
                        else None
                    ),
                    model_path=command.model_path,
                )
                exit_code = self._wait(job, process, command.process_log)
            current = self.repository.get(job.id)
            if current.status == JobStatus.CANCELLING:
                self._fail_if_active(
                    job.id, exit_code, "The job was cancelled.", "cancelled"
                )
                return
            if exit_code != 0:
                self._fail_if_active(
                    job.id, exit_code,
                    f"The worker command exited with code {exit_code}.",
                    "nonzero_exit",
                )
                return
            if job.job_type == "qsar_training":
                self._finish_qsar_training(job, command, exit_code)
                return
            if job.job_type == "qsar_comparison":
                self._finish_qsar_comparison(job, exit_code)
                return
            if job.job_type == "reinforcement_learning":
                self._finish_rl(job, command, exit_code)
                return
            if job.job_type == "transfer_learning":
                self._finish_training(job, command, exit_code)
                return
            self.repository.set_progress(job.id, "validating output")
            created_at = datetime.fromisoformat(job.created_at)
            filename = result_filename(
                comment_slug=job.comment_slug,
                created_at=created_at,
                timezone_name=self.settings.app_timezone,
            )
            final_path = job_dir / "output" / filename
            if command.raw_output is None:
                raise ResultValidationError("The generation output is missing.")
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

    def _finish_training(self, job: JobRecord, command, exit_code: int) -> None:
        output_model = command.output_model_path
        if output_model is None or not output_model.is_file() or output_model.stat().st_size == 0:
            self._fail_if_active(
                job.id,
                exit_code,
                "Transfer Learning completed without a valid model file.",
                "missing_model_output",
            )
            return
        if job.model_id is None or command.config_path is None:
            raise RuntimeError("training model metadata is missing")
        conditions = job.run_conditions
        ModelRepository(self.settings).register(
            model_id=job.model_id,
            name=str(conditions["model_name"]),
            model_type=str(conditions["base_model"]),
            training_data_count=int(conditions["train_count"]),
            validation_data_count=int(conditions["validation_count"]),
            training_conditions=conditions,
            file_path=relative_to_repository(self.settings, output_model),
            config_path=relative_to_repository(self.settings, command.config_path),
            log_path=relative_to_repository(self.settings, command.process_log),
            source_job_id=job.id,
            base_model_key=str(conditions["base_model_key"]),
        )
        self.repository.finish(
            job.id,
            target=JobStatus.COMPLETED,
            exit_code=exit_code,
        )

    def _finish_qsar_training(self, job: JobRecord, command, exit_code: int) -> None:
        conditions = job.run_conditions
        output_model = command.output_model_path
        evaluation_path = confined_path(self.settings.repository_root, conditions["evaluation_path"])
        predictions_path = confined_path(self.settings.repository_root, conditions["predictions_path"])
        dataset_path = confined_path(self.settings.repository_root, conditions["dataset_path"])
        if (
            output_model is None
            or not output_model.is_file()
            or output_model.stat().st_size == 0
            or not evaluation_path.is_file()
            or not predictions_path.is_file()
        ):
            self._fail_if_active(job.id, exit_code, "QSAR training completed without valid artifacts.", "missing_qsar_artifact")
            return
        if job.model_id is None:
            raise RuntimeError("QSAR model metadata is missing")
        evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
        metrics = evaluation["metrics"]
        metadata = evaluation["metadata"]
        QsarModelRepository(self.settings).register(
            model_id=job.model_id,
            name=str(conditions["model_name"]),
            algorithm=str(conditions["algorithm"]),
            task_type=str(conditions["task_type"]),
            target_name=str(conditions["target_name"]),
            target_unit=str(conditions["target_unit"]),
            training_data_count=int(metrics["training_data_count"]),
            test_data_count=int(metrics["test_data_count"]),
            metrics=metrics,
            metadata=metadata,
            file_path=relative_to_repository(self.settings, output_model),
            dataset_path=relative_to_repository(self.settings, dataset_path),
            predictions_path=relative_to_repository(self.settings, predictions_path),
            evaluation_path=relative_to_repository(self.settings, evaluation_path),
            source_job_id=job.id,
            ad_threshold=float(conditions["ad_threshold"]),
        )
        self.repository.set_metrics(job.id, metrics)
        self.repository.finish(job.id, target=JobStatus.COMPLETED, exit_code=exit_code)

    def _finish_qsar_comparison(self, job: JobRecord, exit_code: int) -> None:
        summary_path = confined_path(self.settings.repository_root, job.run_conditions["summary_path"])
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            self._fail_if_active(job.id, exit_code, "QSAR comparison completed without a summary.", "missing_qsar_comparison")
            return
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        self.repository.set_metrics(job.id, summary)
        self.repository.finish(
            job.id,
            target=JobStatus.COMPLETED,
            exit_code=exit_code,
            result_csv_path=relative_to_repository(self.settings, summary_path),
        )

    def _finish_rl(self, job: JobRecord, command, exit_code: int) -> None:
        conditions = job.run_conditions
        output_model = command.output_model_path
        summary_path = confined_path(self.settings.repository_root, conditions["summary_csv_path"])
        if output_model is None or not output_model.is_file() or output_model.stat().st_size == 0:
            self._fail_if_active(job.id, exit_code, "Reinforcement Learning completed without a checkpoint.", "missing_rl_checkpoint")
            return
        if job.model_id is None or command.config_path is None:
            raise RuntimeError("RL model metadata is missing")
        if not summary_path.is_file() or summary_path.stat().st_size == 0:
            self._fail_if_active(job.id, exit_code, "Reinforcement Learning completed without a summary CSV.", "missing_rl_summary")
            return
        audit_path = _write_rl_audit(summary_path)
        ModelRepository(self.settings).register(
            model_id=job.model_id,
            name=str(conditions["model_name"]),
            model_type=str(conditions["model_type"]),
            training_data_count=int(conditions["maximum_steps"]),
            validation_data_count=0,
            training_conditions=conditions,
            file_path=relative_to_repository(self.settings, output_model),
            config_path=relative_to_repository(self.settings, command.config_path),
            log_path=relative_to_repository(self.settings, command.process_log),
            source_job_id=job.id,
            base_model_key=str(conditions["base_model_key"]),
        )
        self.repository.finish(
            job.id,
            target=JobStatus.COMPLETED,
            exit_code=exit_code,
            result_csv_path=relative_to_repository(self.settings, audit_path),
        )

    def _wait(self, job: JobRecord, process: subprocess.Popen, log_path: Path) -> int:
        next_progress_scan = 0.0
        last_epoch = -1
        while True:
            exit_code = process.poll()
            if exit_code is not None:
                return exit_code
            now = time.monotonic()
            if now >= next_progress_scan:
                try:
                    last_epoch = self._update_log_progress(job, log_path, last_epoch)
                except (OSError, RuntimeError, TypeError, ValueError):
                    # Progress extraction is best-effort. A malformed or partial
                    # log line must never turn a successful training process into
                    # a failed job.
                    pass
                next_progress_scan = now + 1.0
            if self.repository.get(job.id).status == JobStatus.CANCELLING:
                self._terminate_process_group(process)
                return process.wait()
            time.sleep(0.05)

    def _update_log_progress(self, job: JobRecord, log_path: Path, last_epoch: int) -> int:
        if job.job_type not in {"transfer_learning", "qsar_training"} or not log_path.is_file():
            return last_epoch
        try:
            with log_path.open("rb") as handle:
                handle.seek(max(0, log_path.stat().st_size - 131_072))
                text = handle.read().decode("utf-8", errors="replace")
        except OSError:
            return last_epoch
        epochs = int(job.run_conditions.get("num_epochs") or job.run_conditions.get("hyperparameters", {}).get("epochs") or 0)
        matches = list(re.finditer(r"Epoch\s+(\d+)\s*/\s*(\d+)", text, re.IGNORECASE))
        if not matches:
            return last_epoch
        raw_epoch = int(matches[-1].group(1))
        raw_total = int(matches[-1].group(2))
        # Chemprop numbers epochs from zero, while the TL runner reports one-based epochs.
        chemprop = job.run_conditions.get("algorithm") == "chemprop_v2"
        epoch = raw_epoch + 1 if chemprop else raw_epoch
        total = epochs or (raw_total + 1 if chemprop else raw_total)
        if epoch <= last_epoch or total < 1:
            return last_epoch
        tail = text[matches[-1].start():]
        number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
        # The trailing boundary prevents checkpoint names such as
        # `val_loss=0.42.ckpt` from being parsed as the invalid float `0.42.`.
        train_matches = re.findall(rf"(?:train_loss_epoch|(?<![\w])loss)\s*[:=]\s*({number})(?![\w.])", tail)
        validation_matches = re.findall(rf"val_loss\s*[:=]\s*({number})(?![\w.])", tail)
        self.repository.set_training_progress(
            job.id,
            epoch=min(epoch, total),
            total_epochs=total,
            training_loss=float(train_matches[-1]) if train_matches else None,
            validation_loss=float(validation_matches[-1]) if validation_matches else None,
        )
        return epoch

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
        if current.job_type in {"transfer_learning", "qsar_training", "reinforcement_learning"}:
            output_relative = str(current.run_conditions.get("output_model_path", ""))
            if output_relative:
                output_model = confined_path(self.settings.repository_root, output_relative)
                if (
                    output_model.is_relative_to(self.settings.model_store_dir)
                    and output_model.is_file()
                    and not output_model.is_symlink()
                ):
                    output_model.unlink()
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
