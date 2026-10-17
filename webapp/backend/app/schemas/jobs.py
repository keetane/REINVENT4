from __future__ import annotations

from datetime import datetime
from typing import Any
from uuid import UUID

from pydantic import BaseModel, Field

from app.models.job import JobStatus


class JobCreateRequest(BaseModel):
    parent_smiles: str = Field(min_length=1, max_length=10_000)
    comment: str = Field(default="", max_length=1_000)
    run_conditions: dict[str, Any] = Field(default_factory=dict)


class JobProgress(BaseModel):
    current: int | None = Field(default=None, ge=0)
    total: int | None = Field(default=None, ge=0)
    percent: float | None = Field(default=None, ge=0, le=100)
    message: str | None = None


class JobResponse(BaseModel):
    id: UUID
    status: JobStatus
    created_at: datetime
    queued_at: datetime
    started_at: datetime | None = None
    finished_at: datetime | None = None
    updated_at: datetime
    comment: str
    parent_smiles: str
    parent_smiles_canonical: str | None = None
    progress: JobProgress = Field(default_factory=JobProgress)
    exit_code: int | None = None
    error_summary: str | None = None


class JobListResponse(BaseModel):
    items: list[JobResponse]
    total: int = Field(ge=0)
    limit: int = Field(ge=1)
    offset: int = Field(ge=0)


class JobLogResponse(BaseModel):
    job_id: UUID
    lines: list[str]
    truncated: bool
