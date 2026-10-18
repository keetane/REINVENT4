from __future__ import annotations

from datetime import datetime
from typing import Any
from uuid import UUID

from pydantic import BaseModel, Field


class ModelResponse(BaseModel):
    id: str
    name: str
    model_type: str
    source: str
    description: str
    created_at: datetime | None = None
    training_data_count: int | None = Field(default=None, ge=0)
    validation_data_count: int | None = Field(default=None, ge=0)
    training_conditions: dict[str, Any] = Field(default_factory=dict)
    base_model_key: str | None = None
    deletable: bool = False


class ModelListResponse(BaseModel):
    items: list[ModelResponse]


class ModelDeleteResponse(BaseModel):
    id: UUID
    deleted: bool
