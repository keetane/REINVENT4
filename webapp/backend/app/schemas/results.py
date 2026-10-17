from __future__ import annotations

from typing import Any
from uuid import UUID

from pydantic import BaseModel, Field


class DescriptorSchema(BaseModel):
    name: str
    source: str
    aliases: list[str] = Field(default_factory=list)
    minimum: float | None = None
    maximum: float | None = None
    missing_count: int = Field(ge=0)


class ResultSchemaResponse(BaseModel):
    job_id: UUID
    columns: list[str]
    descriptors: list[DescriptorSchema]
    total_rows: int = Field(ge=0)
    valid_rows: int = Field(ge=0)
    invalid_smiles: int = Field(ge=0)
    warnings: list[str] = Field(default_factory=list)
    rdkit_version: str


class MoleculeResponse(BaseModel):
    row_number: int = Field(ge=1)
    smiles: str
    canonical_smiles: str
    structure_svg: str
    values: dict[str, Any]


class MoleculeListResponse(BaseModel):
    job_id: UUID
    items: list[MoleculeResponse]
    total: int = Field(ge=0)
    limit: int = Field(ge=1)
    offset: int = Field(ge=0)
    invalid_smiles: int = Field(ge=0)


class NumericSummary(BaseModel):
    count: int = Field(ge=0)
    missing_count: int = Field(ge=0)
    minimum: float | None = None
    maximum: float | None = None
    mean: float | None = None
    median: float | None = None
    standard_deviation: float | None = None
    q1: float | None = None
    q3: float | None = None


class HistogramBin(BaseModel):
    start: float
    end: float
    count: int = Field(ge=0)


class ScatterPoint(BaseModel):
    row_number: int = Field(ge=1)
    x: float
    y: float


class StatisticsResponse(BaseModel):
    job_id: UUID
    total_rows: int = Field(ge=0)
    filtered_rows: int = Field(ge=0)
    invalid_smiles: int = Field(ge=0)
    summaries: dict[str, NumericSummary]
    histograms: dict[str, list[HistogramBin]]
    scatter_x: str | None = None
    scatter_y: str | None = None
    scatter_sampled: bool = False
    scatter: list[ScatterPoint] = Field(default_factory=list)
