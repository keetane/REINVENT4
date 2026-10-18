from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any
from uuid import uuid4
from zoneinfo import ZoneInfo

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.errors import AppError
from app.schemas.qsar import QsarComparisonCreateRequest, QsarDataRequest, QsarJobCreateRequest, QsarRejectedRow
from app.services.csv_data import csv_dict_reader
from app.services.smiles_standardization import (
    StandardizationResult,
    standardization_counts,
    standardization_metadata,
    standardize_smiles,
)
from app.services.workspaces import relative_to_repository


UNIT_TO_MOLAR = {
    "M": 1.0,
    "mM": 1e-3,
    "uM": 1e-6,
    "μM": 1e-6,
    "nM": 1e-9,
    "pM": 1e-12,
}
COMPARATOR = re.compile(r"^\s*[<>]=?")


@dataclass(frozen=True, slots=True)
class ProcessedQsarData:
    columns: list[str]
    total_rows: int
    rows: list[tuple[str, Any]]
    preview: list[dict[str, str]]
    rejected_rows: list[QsarRejectedRow]
    invalid_smiles_count: int
    missing_target_count: int
    invalid_target_count: int
    duplicate_count: int
    conflicting_duplicate_count: int
    standardization_results: list[StandardizationResult]
    standardization_counts: dict[str, int]
    accepted_audit_rows: list[dict[str, Any]]
    original_rows: dict[int, dict[str, str]]


def _target_value(raw: str, request: QsarDataRequest) -> float | str:
    text = str(raw or "").strip()
    if not text:
        raise ValueError("missing_target")
    if COMPARATOR.match(text):
        raise ValueError("censored_target")
    if request.task_type == "classification":
        return text
    try:
        value = float(text)
    except ValueError as error:
        raise ValueError("invalid_target") from error
    if not math.isfinite(value):
        raise ValueError("invalid_target")
    if request.target_transform == "p_activity":
        if value <= 0:
            raise ValueError("nonpositive_concentration")
        value = -math.log10(value * UNIT_TO_MOLAR[request.concentration_unit])
    return value


def process_qsar_data(request: QsarDataRequest, *, max_atoms: int) -> ProcessedQsarData:
    if Path(request.file_name).suffix.lower() != ".csv":
        raise AppError(status_code=422, code="unsupported_qsar_file", message="QSAR学習データはCSV形式で指定してください。")
    reader = csv_dict_reader(request.content)
    columns = list(reader.fieldnames or [])
    if len(columns) < 2:
        raise AppError(status_code=422, code="invalid_qsar_csv", message="CSVにはヘッダーと2列以上のデータが必要です。")
    raw_rows = [dict(row) for row in reader]
    preview = [
        {column: str(row.get(column, "") or "") for column in columns}
        for row in raw_rows[:10]
    ]
    if not request.smiles_column or not request.target_column:
        return ProcessedQsarData(columns, len(raw_rows), [], preview, [], 0, 0, 0, 0, 0, [], standardization_counts([]), [], {})
    if request.smiles_column not in columns or request.target_column not in columns:
        raise AppError(status_code=422, code="qsar_column_not_found", message="選択したSMILES列または目的変数列がCSVにありません。")

    accepted: list[dict[str, Any]] = []
    rejected: list[QsarRejectedRow] = []
    standardization_results: list[StandardizationResult] = []
    invalid_smiles = missing_target = invalid_target = 0
    for row_number, row in enumerate(raw_rows, start=2):
        raw_smiles = str(row.get(request.smiles_column, "") or "").strip()
        raw_target = str(row.get(request.target_column, "") or "").strip()
        standardized = standardize_smiles(raw_smiles)
        standardization_results.append(standardized)
        if not standardized.accepted or not standardized.standardized_smiles:
            invalid_smiles += 1
            rejected.append(QsarRejectedRow(
                row_number=row_number, smiles=raw_smiles, target=raw_target,
                reason=standardized.reason or "invalid_smiles",
                standardized_smiles=standardized.standardized_smiles,
            ))
            continue
        try:
            target = _target_value(raw_target, request)
        except ValueError as error:
            reason = str(error)
            if reason == "missing_target":
                missing_target += 1
            else:
                invalid_target += 1
            rejected.append(QsarRejectedRow(
                row_number=row_number, smiles=raw_smiles, target=raw_target, reason=reason,
                standardized_smiles=standardized.standardized_smiles,
            ))
            continue
        accepted.append({
            "row_number": row_number,
            "smiles": raw_smiles,
            "canonical": standardized.standardized_smiles,
            "molecular_weight": standardized.molecular_weight,
            "formal_charge": standardized.formal_charge,
            "history": standardized.history,
            "target_raw": raw_target,
            "target": target,
        })

    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for item in accepted:
        groups[item["canonical"]].append(item)

    final_rows: list[tuple[str, Any]] = []
    accepted_audit_rows: list[dict[str, Any]] = []
    duplicate_count = 0
    conflicting_count = 0
    for canonical, group in groups.items():
        if len(group) == 1:
            final_rows.append((canonical, group[0]["target"]))
            accepted_audit_rows.append(group[0])
            continue
        duplicate_count += len(group) - 1
        values = [item["target"] for item in group]
        conflict = len(set(values)) > 1
        if conflict:
            conflicting_count += len(group)
            if request.duplicate_strategy == "exclude_conflicts":
                for item in group:
                    rejected.append(QsarRejectedRow(row_number=item["row_number"], smiles=item["smiles"], target=item["target_raw"], reason="conflicting_duplicate", standardized_smiles=canonical))
                continue
            if request.duplicate_strategy == "mean":
                final_rows.append((canonical, sum(float(value) for value in values) / len(values)))
                accepted_audit_rows.append({**group[0], "target": final_rows[-1][1]})
            else:
                final_rows.append((canonical, group[0]["target"]))
                accepted_audit_rows.append(group[0])
        else:
            final_rows.append((canonical, group[0]["target"]))
            accepted_audit_rows.append(group[0])
        for item in group[1:]:
            rejected.append(QsarRejectedRow(row_number=item["row_number"], smiles=item["smiles"], target=item["target_raw"], reason="duplicate", standardized_smiles=canonical))

    return ProcessedQsarData(
        columns=columns,
        total_rows=len(raw_rows),
        rows=final_rows,
        preview=preview,
        rejected_rows=sorted(rejected, key=lambda row: row.row_number),
        invalid_smiles_count=invalid_smiles,
        missing_target_count=missing_target,
        invalid_target_count=invalid_target,
        duplicate_count=duplicate_count,
        conflicting_duplicate_count=conflicting_count,
        standardization_results=standardization_results,
        standardization_counts=standardization_counts(standardization_results),
        accepted_audit_rows=accepted_audit_rows,
        original_rows={index: {column: str(row.get(column, "") or "") for column in columns} for index, row in enumerate(raw_rows, start=2)},
    )


def _write_qsar_audit_artifacts(job_dir: Path, data: ProcessedQsarData) -> tuple[Path, Path]:
    accepted_path = job_dir / "input" / "standardization-accepted.csv"
    with accepted_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["Row", *data.columns, "Standardized_SMILES", "Processed_Target", "MolecularWeight", "FormalCharge", "History"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for item in data.accepted_audit_rows:
            writer.writerow({
                "Row": item["row_number"], **data.original_rows[item["row_number"]],
                "Standardized_SMILES": item["canonical"], "Processed_Target": item["target"],
                "MolecularWeight": item["molecular_weight"], "FormalCharge": item["formal_charge"],
                "History": ";".join(item["history"]),
            })
    rejected_path = job_dir / "input" / "standardization-rejected.csv"
    with rejected_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["Row", *data.columns, "Standardized_SMILES", "Reason"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in data.rejected_rows:
            writer.writerow({
                "Row": row.row_number, **data.original_rows[row.row_number],
                "Standardized_SMILES": row.standardized_smiles or "", "Reason": row.reason,
            })
    return accepted_path, rejected_path


def _next_model_name(settings: Settings, comment: str) -> str:
    date = datetime.now(ZoneInfo(settings.app_timezone)).strftime("%Y%m%d")
    base = f"{date}_{comment}"
    existing = {path.stem for path in settings.qsar_model_store_dir.glob(f"{base}*.*")}
    candidate = base
    serial = 2
    while candidate in existing:
        candidate = f"{base}_{serial}"
        serial += 1
    return candidate


def create_qsar_job(settings: Settings, request: QsarJobCreateRequest) -> JobRecord:
    if request.algorithm == "chemprop_v2" and (not settings.chemprop_runtime_python or not settings.chemprop_runtime_python.is_file()):
        raise AppError(status_code=422, code="chemprop_v2_unavailable", message="Chemprop v2の独立実行環境が見つかりません。")
    if request.algorithm == "qsartuna" and (not settings.qsartuna_runtime_python or not settings.qsartuna_runtime_python.is_file()):
        raise AppError(status_code=422, code="qsartuna_unavailable", message="QSARtunaの独立実行環境が見つかりません。")
    data = process_qsar_data(request, max_atoms=settings.parent_max_atoms)
    minimum = max(6, request.cross_validation_folds + 2)
    if len(data.rows) < minimum:
        raise AppError(status_code=422, code="insufficient_qsar_data", message=f"有効な重複除外済みデータが{minimum}件以上必要です。")
    if request.task_type == "classification":
        counts = Counter(str(target) for _, target in data.rows)
        if len(counts) < 2 or min(counts.values()) < 2:
            raise AppError(status_code=422, code="imbalanced_classification_data", message="分類には2クラス以上かつ各クラス2件以上が必要です。")
        if request.split_method == "random":
            test_count = len(data.rows) - int(len(data.rows) * request.train_fraction)
            train_count = len(data.rows) - test_count
            if min(test_count, train_count) < len(counts):
                raise AppError(
                    status_code=422,
                    code="classification_split_too_small",
                    message="層化分割にはTraining/Testの各側にクラス数以上のデータが必要です。分割比率またはデータ数を調整してください。",
                )

    settings.qsar_model_store_dir.mkdir(parents=True, exist_ok=True)
    model_name = _next_model_name(settings, request.comment)
    model_id = uuid4()
    job_id = uuid4()
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    for name in ("input", "config", "logs", "output"):
        (job_dir / name).mkdir(parents=True, exist_ok=False)

    dataset_path = job_dir / "input" / "training-dataset.csv"
    with dataset_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Canonical_SMILES", "Target"])
        writer.writerows(data.rows)
    accepted_path, rejected_path = _write_qsar_audit_artifacts(job_dir, data)

    extension = ".ckpt" if request.algorithm == "chemprop_v2" else ".pkl" if request.algorithm == "qsartuna" else ".joblib"
    output_model = settings.qsar_model_store_dir / f"{model_name}{extension}"
    output_model.touch(exist_ok=False)
    config_path = job_dir / "config" / "qsar-training.json"
    predictions_path = job_dir / "output" / "predictions.csv"
    evaluation_path = job_dir / "output" / "evaluation.json"
    dataset_hash = hashlib.sha256(dataset_path.read_bytes()).hexdigest()
    target_unit = (
        (f"p{request.activity_type}" if request.activity_type != "custom" else f"p{request.target_column}")
        if request.target_transform == "p_activity"
        else "as_provided"
    )
    conditions = {
        "model_id": str(model_id),
        "model_name": model_name,
        "algorithm": request.algorithm,
        "task_type": request.task_type,
        "target_name": request.target_column,
        "target_unit": target_unit,
        "target_transform": request.target_transform,
        "concentration_unit": request.concentration_unit,
        "activity_type": request.activity_type,
        "duplicate_strategy": request.duplicate_strategy,
        "fingerprint": request.fingerprint,
        "feature_mode": request.feature_mode,
        "morgan_radius": request.morgan_radius,
        "bit_length": request.bit_length,
        "use_chirality": request.use_chirality,
        "split_method": request.split_method,
        "train_fraction": request.train_fraction,
        "cross_validation_folds": request.cross_validation_folds,
        "ad_threshold": request.ad_threshold,
        "device": request.device,
        "advanced": request.advanced,
        "hyperparameters": request.hyperparameters,
        "external_runtime_python": (
            str(settings.chemprop_runtime_python) if request.algorithm == "chemprop_v2"
            else str(settings.qsartuna_runtime_python) if request.algorithm == "qsartuna"
            else None
        ),
        "valid_count": len(data.rows),
        "rejected_count": len(data.rejected_rows),
        "standardization_counts": {**data.standardization_counts, "duplicate_count": data.duplicate_count, "final_accepted_count": len(data.rows)},
        "smiles_standardization": standardization_metadata(),
        "dataset_hash": dataset_hash,
        "source_file_name": Path(request.file_name).name,
        "dataset_path": relative_to_repository(settings, dataset_path),
        "rejected_path": relative_to_repository(settings, rejected_path),
        "standardization_accepted_path": relative_to_repository(settings, accepted_path),
        "standardization_rejected_path": relative_to_repository(settings, rejected_path),
        "output_model_path": relative_to_repository(settings, output_model),
        "predictions_path": relative_to_repository(settings, predictions_path),
        "evaluation_path": relative_to_repository(settings, evaluation_path),
    }
    config_path.write_text(json.dumps(conditions, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    conditions["config_path"] = relative_to_repository(settings, config_path)
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=request.comment,
        parent_smiles="",
        parent_smiles_canonical=None,
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, job_dir),
        job_type="qsar_training",
        model_id=model_id,
    )


def create_qsar_comparison_job(settings: Settings, request: QsarComparisonCreateRequest) -> JobRecord:
    data = process_qsar_data(request, max_atoms=settings.parent_max_atoms)
    minimum = max(6, request.cross_validation_folds + 2)
    if len(data.rows) < minimum:
        raise AppError(status_code=422, code="insufficient_qsar_data", message=f"有効な重複除外済みデータが{minimum}件以上必要です。")
    if request.task_type == "classification":
        counts = Counter(str(target) for _, target in data.rows)
        if len(counts) < 2 or min(counts.values()) < 2:
            raise AppError(status_code=422, code="imbalanced_classification_data", message="分類には2クラス以上かつ各クラス2件以上が必要です。")

    job_id = uuid4()
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    for name in ("input", "config", "logs", "output"):
        (job_dir / name).mkdir(parents=True, exist_ok=False)
    dataset_path = job_dir / "input" / "training-dataset.csv"
    with dataset_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Canonical_SMILES", "Target"])
        writer.writerows(data.rows)
    accepted_path, rejected_path = _write_qsar_audit_artifacts(job_dir, data)

    target_unit = (
        (f"p{request.activity_type}" if request.activity_type != "custom" else f"p{request.target_column}")
        if request.target_transform == "p_activity" else "as_provided"
    )
    common = {
        "task_type": request.task_type,
        "target_name": request.target_column,
        "target_unit": target_unit,
        "target_transform": request.target_transform,
        "concentration_unit": request.concentration_unit,
        "activity_type": request.activity_type,
        "fingerprint": request.fingerprint,
        "feature_mode": request.feature_mode,
        "morgan_radius": request.morgan_radius,
        "bit_length": request.bit_length,
        "use_chirality": request.use_chirality,
        "split_method": request.split_method,
        "train_fraction": request.train_fraction,
        "cross_validation_folds": request.cross_validation_folds,
        "ad_threshold": request.ad_threshold,
        "device": request.device,
        "dataset_path": relative_to_repository(settings, dataset_path),
        "standardization_counts": {**data.standardization_counts, "duplicate_count": data.duplicate_count, "final_accepted_count": len(data.rows)},
        "smiles_standardization": standardization_metadata(),
    }
    models = []
    for algorithm in request.algorithms:
        output_model = job_dir / "output" / f"{algorithm}.joblib"
        predictions_path = job_dir / "output" / f"{algorithm}-predictions.csv"
        evaluation_path = job_dir / "output" / f"{algorithm}-evaluation.json"
        model_config = {
            **common,
            "algorithm": algorithm,
            "model_name": f"{request.comment}_{algorithm}",
            "advanced": request.advanced and algorithm == request.algorithm,
            "hyperparameters": request.hyperparameters if algorithm == request.algorithm else {},
            "optimize_hyperparameters": True,
            "primary_metric": request.primary_metric,
            "output_model_path": relative_to_repository(settings, output_model),
            "predictions_path": relative_to_repository(settings, predictions_path),
            "evaluation_path": relative_to_repository(settings, evaluation_path),
        }
        algorithm_config = job_dir / "config" / f"{algorithm}.json"
        algorithm_config.write_text(json.dumps(model_config, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        models.append({
            "algorithm": algorithm,
            "config_path": relative_to_repository(settings, algorithm_config),
            "output_model_path": model_config["output_model_path"],
            "predictions_path": model_config["predictions_path"],
            "evaluation_path": model_config["evaluation_path"],
        })
    summary_path = job_dir / "output" / "comparison.json"
    comparison_config = job_dir / "config" / "qsar-comparison.json"
    comparison_config.write_text(json.dumps({"primary_metric": request.primary_metric, "summary_path": relative_to_repository(settings, summary_path), "models": models}, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    conditions = {
        **common,
        "comment": request.comment,
        "algorithms": request.algorithms,
        "primary_metric": request.primary_metric,
        "models": models,
        "valid_count": len(data.rows),
        "rejected_count": len(data.rejected_rows),
        "dataset_hash": hashlib.sha256(dataset_path.read_bytes()).hexdigest(),
        "dataset_path": relative_to_repository(settings, dataset_path),
        "rejected_path": relative_to_repository(settings, rejected_path),
        "standardization_accepted_path": relative_to_repository(settings, accepted_path),
        "standardization_rejected_path": relative_to_repository(settings, rejected_path),
        "summary_path": relative_to_repository(settings, summary_path),
        "config_path": relative_to_repository(settings, comparison_config),
    }
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=request.comment,
        parent_smiles="",
        parent_smiles_canonical=None,
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, job_dir),
        job_type="qsar_comparison",
    )
