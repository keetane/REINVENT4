from __future__ import annotations

from typing import Literal
from uuid import UUID

from pydantic import BaseModel, Field, model_validator


class RlQsarComponent(BaseModel):
    model_id: UUID
    weight: float = Field(default=1.0, gt=0, le=100)
    transform: Literal["sigmoid", "reverse_sigmoid", "double_sigmoid", "none"] = "sigmoid"
    low: float = 0.0
    high: float = 1.0
    k: float = Field(default=0.5, gt=0, le=100)
    class_index: int = Field(default=1, ge=0, le=100)
    ad_behavior: Literal["penalize", "allow"] = "penalize"

    @model_validator(mode="after")
    def check_bounds(self):
        if self.transform != "none" and self.low >= self.high:
            raise ValueError("score transform low must be lower than high")
        return self


class RlPropertyComponent(BaseModel):
    component: Literal[
        "QED", "MolecularWeight", "SlogP", "TPSA", "HBondDonors",
        "HBondAcceptors", "NumRotBond", "SAScore", "TanimotoSimilarity",
        "CustomAlerts",
    ]
    weight: float = Field(default=1.0, gt=0, le=100)
    transform: Literal["sigmoid", "reverse_sigmoid", "double_sigmoid", "none"] = "none"
    low: float = 0.0
    high: float = 1.0
    k: float = Field(default=0.5, gt=0, le=100)
    reference_smiles: list[str] = Field(default_factory=list, max_length=100)
    smarts: list[str] = Field(default_factory=list, max_length=100)

    @model_validator(mode="after")
    def check_component(self):
        if self.transform != "none" and self.low >= self.high:
            raise ValueError("score transform low must be lower than high")
        if self.component == "TanimotoSimilarity" and not self.reference_smiles:
            raise ValueError("TanimotoSimilarity requires reference_smiles")
        if self.component == "CustomAlerts" and not self.smarts:
            raise ValueError("CustomAlerts requires SMARTS patterns")
        return self


class RlJobCreateRequest(BaseModel):
    generation_model_id: str = Field(min_length=1, max_length=100)
    comment: str = Field(pattern=r"^[A-Za-z0-9_-]{1,50}$")
    device: Literal["cpu", "cuda", "mps"] = "mps"
    input_smiles_content: str = Field(default="", max_length=10_000_000)
    qsar_components: list[RlQsarComponent] = Field(default_factory=list, max_length=20)
    property_components: list[RlPropertyComponent] = Field(default_factory=list, max_length=20)
    aggregation: Literal["geometric_mean", "arithmetic_mean"] = "geometric_mean"
    batch_size: int = Field(default=64, ge=1, le=4096)
    sigma: int = Field(default=128, ge=1, le=10000)
    learning_rate: float = Field(default=0.0001, ge=1e-8, le=1.0)
    minimum_steps: int = Field(default=25, ge=0, le=100000)
    maximum_steps: int = Field(default=300, ge=1, le=100000)
    target_mean_score: float = Field(default=0.8, ge=0, le=1)
    diversity_filter: bool = True
    diversity_bucket_size: int = Field(default=25, ge=1, le=10000)
    diversity_min_score: float = Field(default=0.4, ge=0, le=1)

    @model_validator(mode="after")
    def check_job(self):
        if self.minimum_steps > self.maximum_steps:
            raise ValueError("minimum_steps must not exceed maximum_steps")
        if not self.qsar_components and not self.property_components:
            raise ValueError("at least one scoring component is required")
        return self


class RlProgressRequest(BaseModel):
    step: int = Field(ge=0)
    stage: int = Field(ge=1)
    components: dict[str, float | None] = Field(default_factory=dict)
    learning: dict[str, float | None] = Field(default_factory=dict)
    fraction_valid_smiles: float | None = Field(default=None, ge=0, le=1)
    collected_smiles_in_memory: int | None = Field(default=None, ge=0, alias="collected smiles in memory")


class QsarPredictionRequest(BaseModel):
    model_id: UUID
    smiles: list[str] = Field(min_length=1, max_length=1000)


class QsarPredictionResponse(BaseModel):
    model_id: UUID
    predictions: list[dict]


class RlPreviewRequest(BaseModel):
    smiles: list[str] = Field(min_length=1, max_length=1000)
    qsar_components: list[RlQsarComponent] = Field(default_factory=list, max_length=20)
    property_components: list[RlPropertyComponent] = Field(default_factory=list, max_length=20)
    aggregation: Literal["geometric_mean", "arithmetic_mean"] = "geometric_mean"

    @model_validator(mode="after")
    def check_components(self):
        if not self.qsar_components and not self.property_components:
            raise ValueError("at least one scoring component is required")
        return self


class RlPreviewResponse(BaseModel):
    rows: list[dict]
