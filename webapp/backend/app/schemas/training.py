from __future__ import annotations

from pydantic import BaseModel, Field, field_validator, model_validator

from app.constants import MAX_DATA_UPLOAD_BYTES


class TrainingDataRequest(BaseModel):
    file_name: str = Field(min_length=1, max_length=255)
    content: str = Field(min_length=1, max_length=MAX_DATA_UPLOAD_BYTES)
    smiles_column: str | None = Field(default=None, max_length=255)
    base_model: str | None = None

    @field_validator("content")
    @classmethod
    def validate_content_bytes(cls, value: str) -> str:
        if len(value.encode("utf-8")) > MAX_DATA_UPLOAD_BYTES:
            raise ValueError("uploaded data must not exceed 50 MiB")
        return value


class TrainingDataInspectionResponse(BaseModel):
    columns: list[str]
    total_rows: int = Field(ge=0)
    valid_count: int = Field(ge=0)
    invalid_count: int = Field(ge=0)
    duplicate_count: int = Field(ge=0)
    standardized_count: int = Field(ge=0)
    fragment_extracted_count: int = Field(ge=0)
    neutralized_count: int = Field(ge=0)
    charged_structure_retained_count: int = Field(ge=0)
    low_molecular_weight_excluded_count: int = Field(ge=0)
    model_compatibility_error_count: int = Field(ge=0)
    standardization_version: str


class TrainingJobCreateRequest(TrainingDataRequest):
    base_model: str
    comment: str = Field(pattern=r"^[A-Za-z0-9_-]{1,50}$")
    train_fraction: float = Field(default=0.8, ge=0.5, le=0.95)
    num_epochs: int = Field(default=3, ge=1, le=1000)
    batch_size: int = Field(default=50, ge=1, le=4096)
    learning_rate: float = Field(default=0.0001, gt=0, le=1)
    tanimoto_lower: float = Field(default=0.7, ge=0, le=1)
    tanimoto_upper: float = Field(default=1.0, ge=0, le=1)
    device: str = "mps"

    @model_validator(mode="after")
    def validate_options(self):
        if self.base_model not in {"reinvent", "mol2mol"}:
            raise ValueError("base_model must be reinvent or mol2mol")
        if self.device not in {"cpu", "cuda", "mps"}:
            raise ValueError("device must be cpu, cuda, or mps")
        if self.tanimoto_lower > self.tanimoto_upper:
            raise ValueError("tanimoto_lower must not exceed tanimoto_upper")
        return self


class TrainingProgressRequest(BaseModel):
    epoch: int = Field(ge=1)
    learning_mean: dict[str, float | None]
