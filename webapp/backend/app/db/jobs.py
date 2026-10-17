from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any
from uuid import UUID

from app.config import Settings
from app.db.connection import connect
from app.models.job import JobStatus, require_transition


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True, slots=True)
class JobRecord:
    id: UUID
    status: JobStatus
    created_at: str
    queued_at: str
    started_at: str | None
    finished_at: str | None
    updated_at: str
    comment_original: str
    comment_slug: str
    parent_smiles_input: str
    parent_smiles_canonical: str | None
    run_conditions: dict[str, Any]
    worker_id: str | None
    process_id: int | None
    exit_code: int | None
    progress_current: int | None
    progress_total: int | None
    progress_percent: float | None
    progress_message: str | None
    job_dir: str
    log_path: str | None
    result_csv_path: str | None
    result_sdf_path: str | None
    error_summary: str | None
    error_kind: str | None
    cancel_requested_at: str | None
    config_snapshot_path: str | None
    model_path: str | None
    command_metadata: dict[str, Any] | None


def _record(row: sqlite3.Row) -> JobRecord:
    return JobRecord(
        id=UUID(row["id"]),
        status=JobStatus(row["status"]),
        created_at=row["created_at"],
        queued_at=row["queued_at"],
        started_at=row["started_at"],
        finished_at=row["finished_at"],
        updated_at=row["updated_at"],
        comment_original=row["comment_original"],
        comment_slug=row["comment_slug"],
        parent_smiles_input=row["parent_smiles_input"],
        parent_smiles_canonical=row["parent_smiles_canonical"],
        run_conditions=json.loads(row["run_conditions_json"]),
        worker_id=row["worker_id"],
        process_id=row["process_id"],
        exit_code=row["exit_code"],
        progress_current=row["progress_current"],
        progress_total=row["progress_total"],
        progress_percent=row["progress_percent"],
        progress_message=row["progress_message"],
        job_dir=row["job_dir"],
        log_path=row["log_path"],
        result_csv_path=row["result_csv_path"],
        result_sdf_path=row["result_sdf_path"],
        error_summary=row["error_summary"],
        error_kind=row["error_kind"],
        cancel_requested_at=row["cancel_requested_at"],
        config_snapshot_path=row["config_snapshot_path"],
        model_path=row["model_path"],
        command_metadata=(
            json.loads(row["command_metadata_json"])
            if row["command_metadata_json"]
            else None
        ),
    )


class JobRepository:
    def __init__(self, settings: Settings):
        self.settings = settings

    def create(
        self,
        *,
        job_id: UUID,
        comment_original: str,
        comment_slug: str,
        parent_smiles: str,
        parent_smiles_canonical: str | None = None,
        run_conditions: dict[str, Any],
        job_dir: str,
    ) -> JobRecord:
        now = utc_now()
        with connect(self.settings) as connection:
            connection.execute(
                """
                INSERT INTO jobs (
                    id, status, created_at, queued_at, updated_at,
                    comment_original, comment_slug, parent_smiles_input,
                    parent_smiles_canonical, run_conditions_json, job_dir,
                    progress_message
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    str(job_id),
                    JobStatus.QUEUED.value,
                    now,
                    now,
                    now,
                    comment_original,
                    comment_slug,
                    parent_smiles,
                    parent_smiles_canonical,
                    json.dumps(run_conditions, ensure_ascii=False, sort_keys=True),
                    job_dir,
                    "queued",
                ),
            )
        return self.get(job_id)

    def get(self, job_id: UUID) -> JobRecord:
        with connect(self.settings) as connection:
            row = connection.execute(
                "SELECT * FROM jobs WHERE id = ?", (str(job_id),)
            ).fetchone()
        if row is None:
            raise KeyError(str(job_id))
        return _record(row)

    def list(self, *, limit: int, offset: int) -> tuple[list[JobRecord], int]:
        with connect(self.settings) as connection:
            total = connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
            rows = connection.execute(
                "SELECT * FROM jobs ORDER BY created_at DESC LIMIT ? OFFSET ?",
                (limit, offset),
            ).fetchall()
        return [_record(row) for row in rows], total

    def claim_next(self, worker_id: str) -> JobRecord | None:
        now = utc_now()
        with connect(self.settings) as connection:
            connection.execute("BEGIN IMMEDIATE")
            try:
                row = connection.execute(
                    """
                    SELECT id FROM jobs
                    WHERE status = ?
                    ORDER BY queued_at, id
                    LIMIT 1
                    """,
                    (JobStatus.QUEUED.value,),
                ).fetchone()
                if row is None:
                    connection.execute("COMMIT")
                    return None
                updated = connection.execute(
                    """
                    UPDATE jobs
                    SET status = ?, worker_id = ?, started_at = ?, updated_at = ?,
                        progress_message = ?
                    WHERE id = ? AND status = ?
                    """,
                    (
                        JobStatus.RUNNING.value,
                        worker_id,
                        now,
                        now,
                        "preparing",
                        row["id"],
                        JobStatus.QUEUED.value,
                    ),
                ).rowcount
                connection.execute("COMMIT")
            except Exception:
                connection.execute("ROLLBACK")
                raise
        if updated != 1:
            return None
        return self.get(UUID(row["id"]))

    def request_cancel(self, job_id: UUID) -> JobRecord:
        current = self.get(job_id)
        now = utc_now()
        if current.status == JobStatus.QUEUED:
            target = JobStatus.CANCELLED
            finished_at = now
        elif current.status == JobStatus.RUNNING:
            target = JobStatus.CANCELLING
            finished_at = None
        elif current.status in {JobStatus.CANCELLING, JobStatus.CANCELLED}:
            return current
        else:
            require_transition(current.status, JobStatus.CANCELLING)
            raise AssertionError("unreachable")
        require_transition(current.status, target)
        with connect(self.settings) as connection:
            updated = connection.execute(
                """
                UPDATE jobs
                SET status = ?, cancel_requested_at = ?, finished_at = ?,
                    updated_at = ?, progress_message = ?
                WHERE id = ? AND status = ?
                """,
                (
                    target.value,
                    now,
                    finished_at,
                    now,
                    target.value,
                    str(job_id),
                    current.status.value,
                ),
            ).rowcount
        if updated != 1:
            return self.get(job_id)
        return self.get(job_id)

    def set_process(
        self,
        job_id: UUID,
        *,
        process_id: int,
        log_path: str,
        command: list[str],
        config_path: str | None = None,
        model_path: str | None = None,
    ) -> None:
        with connect(self.settings) as connection:
            connection.execute(
                """
                UPDATE jobs
                SET process_id = ?, log_path = ?, command_metadata_json = ?,
                    config_snapshot_path = ?, model_path = ?, updated_at = ?,
                    progress_message = ?
                WHERE id = ? AND status = ?
                """,
                (
                    process_id,
                    log_path,
                    json.dumps({"argv": command}, ensure_ascii=False),
                    config_path,
                    model_path,
                    utc_now(),
                    "running",
                    str(job_id),
                    JobStatus.RUNNING.value,
                ),
            )

    def set_progress(self, job_id: UUID, message: str) -> None:
        with connect(self.settings) as connection:
            connection.execute(
                """
                UPDATE jobs SET progress_message = ?, updated_at = ?
                WHERE id = ? AND status IN (?, ?)
                """,
                (
                    message,
                    utc_now(),
                    str(job_id),
                    JobStatus.RUNNING.value,
                    JobStatus.CANCELLING.value,
                ),
            )

    def finish(
        self,
        job_id: UUID,
        *,
        target: JobStatus,
        exit_code: int | None,
        result_csv_path: str | None = None,
        error_summary: str | None = None,
        error_kind: str | None = None,
    ) -> JobRecord:
        current = self.get(job_id)
        require_transition(current.status, target)
        now = utc_now()
        with connect(self.settings) as connection:
            updated = connection.execute(
                """
                UPDATE jobs
                SET status = ?, finished_at = ?, updated_at = ?, exit_code = ?,
                    result_csv_path = ?, error_summary = ?, error_kind = ?,
                    progress_message = ?, progress_percent = ?
                WHERE id = ? AND status = ?
                """,
                (
                    target.value,
                    now,
                    now,
                    exit_code,
                    result_csv_path,
                    error_summary,
                    error_kind,
                    target.value,
                    100.0 if target == JobStatus.COMPLETED else None,
                    str(job_id),
                    current.status.value,
                ),
            ).rowcount
        if updated != 1:
            raise RuntimeError("job changed while finishing")
        return self.get(job_id)

    def interrupt(self, job_id: UUID, message: str) -> JobRecord:
        current = self.get(job_id)
        return self.finish(
            job_id,
            target=JobStatus.INTERRUPTED,
            exit_code=current.exit_code,
            error_summary=message,
            error_kind="worker_interrupted",
        )

    def active_jobs(self) -> list[JobRecord]:
        with connect(self.settings) as connection:
            rows = connection.execute(
                "SELECT * FROM jobs WHERE status IN (?, ?)",
                (JobStatus.RUNNING.value, JobStatus.CANCELLING.value),
            ).fetchall()
        return [_record(row) for row in rows]
