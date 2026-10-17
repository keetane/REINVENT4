from __future__ import annotations

from collections import deque
from pathlib import Path
from uuid import UUID

from fastapi import APIRouter, Query, Request, status

from app.db.jobs import JobRecord, JobRepository
from app.errors import AppError
from app.models.job import TransitionError
from app.schemas.jobs import (
    JobCreateRequest,
    JobListResponse,
    JobLogResponse,
    JobProgress,
    JobResponse,
)
from app.services.jobs import create_job
from app.services.naming import confined_path
from app.services.log_sanitizer import sanitize_log_lines


router = APIRouter(prefix="/jobs", tags=["jobs"])
MAX_LOG_BYTES = 128 * 1024


def _response(job: JobRecord) -> JobResponse:
    return JobResponse(
        id=job.id,
        status=job.status,
        created_at=job.created_at,
        queued_at=job.queued_at,
        started_at=job.started_at,
        finished_at=job.finished_at,
        updated_at=job.updated_at,
        comment=job.comment_original,
        parent_smiles=job.parent_smiles_input,
        parent_smiles_canonical=job.parent_smiles_canonical,
        progress=JobProgress(
            current=job.progress_current,
            total=job.progress_total,
            percent=job.progress_percent,
            message=job.progress_message,
        ),
        exit_code=job.exit_code,
        error_summary=job.error_summary,
    )


def _get(repository: JobRepository, job_id: UUID) -> JobRecord:
    try:
        return repository.get(job_id)
    except KeyError as error:
        raise AppError(
            status_code=404,
            code="job_not_found",
            message="The requested job does not exist.",
        ) from error


@router.post("", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED)
def post_job(request: Request, payload: JobCreateRequest) -> JobResponse:
    return _response(create_job(request.app.state.settings, payload))


@router.get("", response_model=JobListResponse)
def get_jobs(
    request: Request,
    limit: int = Query(default=25, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
) -> JobListResponse:
    jobs, total = JobRepository(request.app.state.settings).list(limit=limit, offset=offset)
    return JobListResponse(
        items=[_response(job) for job in jobs], total=total, limit=limit, offset=offset
    )


@router.get("/{job_id}", response_model=JobResponse)
def get_job(request: Request, job_id: UUID) -> JobResponse:
    return _response(_get(JobRepository(request.app.state.settings), job_id))


@router.post("/{job_id}/cancel", response_model=JobResponse)
def cancel_job(request: Request, job_id: UUID) -> JobResponse:
    repository = JobRepository(request.app.state.settings)
    _get(repository, job_id)
    try:
        return _response(repository.request_cancel(job_id))
    except TransitionError as error:
        raise AppError(
            status_code=409,
            code="job_not_cancellable",
            message="The job can no longer be cancelled.",
        ) from error


@router.get("/{job_id}/logs", response_model=JobLogResponse)
def get_job_logs(
    request: Request,
    job_id: UUID,
    tail: int = Query(default=100, ge=1, le=1_000),
) -> JobLogResponse:
    settings = request.app.state.settings
    job = _get(JobRepository(settings), job_id)
    if not job.log_path:
        return JobLogResponse(job_id=job_id, lines=[], truncated=False)
    try:
        log_path = confined_path(settings.repository_root, job.log_path)
    except ValueError as error:
        raise AppError(
            status_code=500,
            code="invalid_log_path",
            message="The job log is unavailable.",
        ) from error
    if not log_path.is_file() or log_path.is_symlink():
        return JobLogResponse(job_id=job_id, lines=[], truncated=False)
    data = _bounded_tail(log_path, MAX_LOG_BYTES)
    decoded = data.decode("utf-8", errors="replace")
    all_lines = decoded.splitlines()
    visible_lines = list(deque(all_lines, maxlen=tail))
    return JobLogResponse(
        job_id=job_id,
        lines=sanitize_log_lines(visible_lines, settings),
        truncated=log_path.stat().st_size > len(data) or len(all_lines) > tail,
    )


def _bounded_tail(path: Path, max_bytes: int) -> bytes:
    with path.open("rb") as handle:
        handle.seek(0, 2)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        return handle.read(max_bytes)
