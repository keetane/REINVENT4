from __future__ import annotations

from datetime import datetime
from typing import Any
from uuid import UUID

from pydantic import BaseModel, Field, model_validator

from app.models.job import JobStatus


class JobCreateRequest(BaseModel):
    parent_smiles: str = Field(default="", max_length=10_000)
    comment: str = Field(default="", max_length=1_000)
    run_conditions: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def require_parent_unless_model_can_resolve_it(self):
        model_id = str(self.run_conditions.get("model_id", ""))
        if not self.parent_smiles.strip() and model_id not in {"standard:reinvent"}:
            if not model_id or model_id.startswith("standard:"):
                raise ValueError("parent_smiles is required for this model")
        return self


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
    job_type: str = "generation"
    model_id: UUID | None = None
    metrics: dict[str, Any] = Field(default_factory=dict)


class JobListResponse(BaseModel):
    items: list[JobResponse]
    total: int = Field(ge=0)
    limit: int = Field(ge=1)
    offset: int = Field(ge=0)


class JobLogResponse(BaseModel):
    job_id: UUID
    lines: list[str]
    truncated: bool
