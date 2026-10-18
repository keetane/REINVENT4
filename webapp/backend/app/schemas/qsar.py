from __future__ import annotations

from datetime import datetime
from typing import Any, Literal
from uuid import UUID

from pydantic import BaseModel, Field, field_validator, model_validator

from app.constants import MAX_DATA_UPLOAD_BYTES


TaskType = Literal["regression", "classification"]
Algorithm = Literal[
    "random_forest", "svm", "gradient_boosting", "extra_trees",
    "lightgbm", "xgboost", "catboost", "kernel_ridge", "gaussian_process",
    "ridge", "elastic_net", "logistic_regression", "mlp",
    "chemprop_v2", "qsartuna",
]


class QsarCapability(BaseModel):
    key: str
    label: str
    available: bool
    reason: str | None = None
    version: str | None = None
    task_types: list[TaskType] = Field(default_factory=lambda: ["regression", "classification"])


class QsarCapabilitiesResponse(BaseModel):
    algorithms: list[QsarCapability]
    devices: list[QsarCapability]


class QsarDataRequest(BaseModel):
    file_name: str = Field(min_length=1, max_length=255)
    content: str = Field(min_length=1, max_length=MAX_DATA_UPLOAD_BYTES)
    smiles_column: str | None = Field(default=None, max_length=255)
    target_column: str | None = Field(default=None, max_length=255)
    task_type: TaskType = "regression"
    target_transform: Literal["none", "p_activity"] = "none"
    concentration_unit: Literal["M", "mM", "uM", "μM", "nM", "pM"] = "nM"
    activity_type: Literal["IC50", "EC50", "Ki", "Kd", "custom"] = "IC50"
    duplicate_strategy: Literal["exclude_conflicts", "keep_first", "mean"] = "exclude_conflicts"

    @field_validator("content")
    @classmethod
    def validate_content_bytes(cls, value: str) -> str:
        if len(value.encode("utf-8")) > MAX_DATA_UPLOAD_BYTES:
            raise ValueError("uploaded data must not exceed 50 MiB")
        return value

    @model_validator(mode="after")
    def validate_transform(self):
        if self.task_type == "classification" and self.target_transform != "none":
            raise ValueError("classification targets cannot use concentration conversion")
        if self.task_type == "classification" and self.duplicate_strategy == "mean":
            raise ValueError("classification targets cannot average duplicate values")
        return self


class QsarRejectedRow(BaseModel):
    row_number: int
    smiles: str
    target: str
    reason: str
    standardized_smiles: str | None = None


class QsarDataInspectionResponse(BaseModel):
    columns: list[str]
    total_rows: int = Field(ge=0)
    valid_count: int = Field(ge=0)
    invalid_smiles_count: int = Field(ge=0)
    missing_target_count: int = Field(ge=0)
    invalid_target_count: int = Field(ge=0)
    duplicate_count: int = Field(ge=0)
    conflicting_duplicate_count: int = Field(ge=0)
    standardized_count: int = Field(ge=0)
    fragment_extracted_count: int = Field(ge=0)
    neutralized_count: int = Field(ge=0)
    charged_structure_retained_count: int = Field(ge=0)
    low_molecular_weight_excluded_count: int = Field(ge=0)
    final_accepted_count: int = Field(ge=0)
    standardization_version: str
    preview: list[dict[str, str]]
    rejected_rows: list[QsarRejectedRow]


class QsarJobCreateRequest(QsarDataRequest):
    algorithm: Algorithm = "random_forest"
    comment: str = Field(pattern=r"^[A-Za-z0-9_-]{1,50}$")
    fingerprint: Literal["morgan", "rdkit", "maccs"] = "morgan"
    feature_mode: Literal["fingerprint", "descriptors", "both"] = "fingerprint"
    morgan_radius: int = Field(default=2, ge=1, le=5)
    bit_length: int = Field(default=2048, ge=128, le=8192)
    use_chirality: bool = True
    split_method: Literal["scaffold", "random"] = "scaffold"
    train_fraction: float = Field(default=0.8, ge=0.5, le=0.95)
    cross_validation_folds: int = Field(default=3, ge=2, le=10)
    ad_threshold: float = Field(default=0.4, ge=0, le=1)
    device: Literal["cpu", "cuda", "mps"] = "cpu"
    advanced: bool = False
    hyperparameters: dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="after")
    def validate_algorithm(self):
        allowed = {
            "random_forest": {"n_estimators", "max_depth", "min_samples_split", "min_samples_leaf"},
            "svm": {"C", "kernel", "gamma", "epsilon"},
            "gradient_boosting": {"n_estimators", "learning_rate", "max_depth"},
            "extra_trees": {"n_estimators", "max_depth", "min_samples_leaf"},
            "lightgbm": {"n_estimators", "learning_rate", "max_depth", "num_leaves", "min_child_samples", "subsample", "colsample_bytree"},
            "xgboost": {"n_estimators", "learning_rate", "max_depth", "min_child_weight", "subsample", "colsample_bytree", "reg_alpha", "reg_lambda"},
            "catboost": {"iterations", "learning_rate", "depth", "l2_leaf_reg"},
            "kernel_ridge": {"alpha", "kernel", "gamma", "degree"},
            "gaussian_process": {"kernel"},
            "ridge": {"alpha"},
            "elastic_net": {"alpha", "l1_ratio", "max_iter"},
            "logistic_regression": {"C", "max_iter"},
            "mlp": {"hidden_layer_sizes", "alpha", "learning_rate_init", "max_iter"},
            "chemprop_v2": {"epochs", "batch_size", "learning_rate", "patience"},
            "qsartuna": {"n_trials", "n_startup_trials", "timeout_seconds"},
        }[self.algorithm]
        if set(self.hyperparameters) - allowed:
            raise ValueError("unsupported hyperparameter for the selected algorithm")
        if self.algorithm == "svm" and self.hyperparameters.get("kernel", "rbf") not in {"linear", "poly", "rbf", "sigmoid"}:
            raise ValueError("unsupported SVM kernel")
        if self.algorithm == "kernel_ridge" and self.hyperparameters.get("kernel", "rbf") not in {"linear", "poly", "rbf", "laplacian", "sigmoid"}:
            raise ValueError("unsupported Kernel Ridge kernel")
        if self.algorithm == "gaussian_process" and self.hyperparameters.get("kernel", "rbf") not in {"rbf", "matern"}:
            raise ValueError("unsupported Gaussian Process kernel")
        regression_only = {"kernel_ridge", "ridge", "elastic_net"}
        if self.task_type == "classification" and self.algorithm in regression_only:
            raise ValueError("the selected algorithm supports regression only")
        if self.task_type == "regression" and self.algorithm == "logistic_regression":
            raise ValueError("Logistic Regression supports classification only")
        numeric_bounds = {
            "n_estimators": (1, 10000), "max_depth": (-1, 1000),
            "min_samples_split": (2, 1000000), "min_samples_leaf": (1, 1000000),
            "C": (1e-12, 1e12), "epsilon": (0, 1e12),
            "learning_rate": (1e-12, 1.0),
            "iterations": (1, 100000), "depth": (1, 16), "num_leaves": (2, 100000),
            "min_child_samples": (1, 1000000), "min_child_weight": (0, 1e12),
            "subsample": (1e-12, 1.0), "colsample_bytree": (1e-12, 1.0),
            "reg_alpha": (0, 1e12), "reg_lambda": (0, 1e12), "l2_leaf_reg": (0, 1e12),
            "alpha": (0, 1e12), "l1_ratio": (0, 1.0), "max_iter": (1, 1000000),
            "learning_rate_init": (1e-12, 1.0), "degree": (1, 100), "gamma": (1e-12, 1e12),
            "epochs": (1, 10000), "batch_size": (1, 100000), "patience": (1, 1000),
            "n_trials": (1, 100000), "n_startup_trials": (0, 100000), "timeout_seconds": (1, 604800),
        }
        for key, value in self.hyperparameters.items():
            if key in {"kernel", "gamma"} and isinstance(value, str):
                if key == "gamma" and value not in {"scale", "auto"}:
                    raise ValueError("unsupported SVM gamma")
                continue
            if key == "hidden_layer_sizes":
                if not isinstance(value, list) or not value or any(isinstance(item, bool) or not 1 <= int(item) <= 10000 for item in value):
                    raise ValueError("hidden_layer_sizes must be a non-empty list of positive integers")
                continue
            if key in numeric_bounds:
                low, high = numeric_bounds[key]
                if isinstance(value, bool) or not low <= float(value) <= high:
                    raise ValueError(f"{key} is outside its supported range")
                if key == "max_depth" and float(value) == -1 and self.algorithm != "lightgbm":
                    raise ValueError("max_depth=-1 is supported only by LightGBM")
        return self


class QsarComparisonCreateRequest(QsarJobCreateRequest):
    algorithms: list[Algorithm] = Field(min_length=2, max_length=15)
    primary_metric: Literal["rmse", "mae", "r2", "roc_auc", "pr_auc"] = "rmse"

    @model_validator(mode="after")
    def validate_comparison(self):
        if len(set(self.algorithms)) != len(self.algorithms):
            raise ValueError("comparison algorithms must be unique")
        unavailable = {"chemprop_v2", "qsartuna"}
        if set(self.algorithms) & unavailable:
            raise ValueError("the selected comparison includes an unavailable algorithm")
        regression_only = {"kernel_ridge", "ridge", "elastic_net"}
        if self.task_type == "classification" and set(self.algorithms) & regression_only:
            raise ValueError("the comparison includes a regression-only algorithm")
        if self.task_type == "regression" and "logistic_regression" in self.algorithms:
            raise ValueError("the comparison includes a classification-only algorithm")
        valid_metrics = {"rmse", "mae", "r2"} if self.task_type == "regression" else {"roc_auc", "pr_auc"}
        if self.primary_metric not in valid_metrics:
            raise ValueError("primary metric does not match the task type")
        return self


class QsarComparisonRegistrationRequest(BaseModel):
    algorithm: Algorithm
    comment: str = Field(pattern=r"^[A-Za-z0-9_-]{1,50}$")


class QsarModelResponse(BaseModel):
    id: UUID
    name: str
    algorithm: str
    task_type: str
    target_name: str
    target_unit: str
    created_at: datetime
    training_data_count: int
    test_data_count: int
    metrics: dict[str, Any]
    metadata: dict[str, Any]
    ad_threshold: float
    in_use: bool = False


class QsarModelListResponse(BaseModel):
    items: list[QsarModelResponse]


class QsarModelDeleteResponse(BaseModel):
    id: UUID
    deleted: bool
