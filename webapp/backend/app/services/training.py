from __future__ import annotations

import json
import random
import csv
import os
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from uuid import uuid4
from zoneinfo import ZoneInfo

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.errors import AppError
from app.schemas.training import TrainingDataRequest, TrainingJobCreateRequest
from app.services.csv_data import csv_dict_reader
from app.services.smiles_standardization import (
    StandardizationResult,
    standardization_counts,
    standardization_metadata,
    standardize_smiles,
)
from app.services.workspaces import relative_to_repository


@dataclass(frozen=True, slots=True)
class TrainingDataInspection:
    columns: list[str]
    total_rows: int
    valid_smiles: list[str]
    invalid_count: int
    duplicate_count: int
    results: list[StandardizationResult]
    counts: dict[str, int]
    model_compatibility_error_count: int = 0
    compatibility_errors: tuple[str, ...] = ()
    compatibility_rejections: tuple[tuple[str, str], ...] = ()


def _extract_values(request: TrainingDataRequest) -> tuple[list[str], list[str]]:
    suffix = Path(request.file_name).suffix.lower()
    if suffix == ".csv":
        reader = csv_dict_reader(request.content)
        columns = list(reader.fieldnames or [])
        if not columns:
            raise AppError(
                status_code=422,
                code="invalid_training_csv",
                message="CSVのヘッダーを読み取れませんでした。",
            )
        if not request.smiles_column:
            return columns, []
        if request.smiles_column not in columns:
            raise AppError(
                status_code=422,
                code="smiles_column_not_found",
                message="選択したSMILES列がCSVにありません。",
            )
        return columns, [row.get(request.smiles_column, "") for row in reader]

    if suffix not in {".smi", ".smiles", ".txt"}:
        raise AppError(
            status_code=422,
            code="unsupported_training_file",
            message="SMILES、SMI、TXT、CSVファイルを選択してください。",
        )
    values = []
    for line in request.content.splitlines():
        stripped = line.strip()
        if stripped:
            values.append(stripped.split()[0])
    return [], values


def _model_compatibility(
    settings: Settings | None, base_model: str | None, smiles: list[str]
) -> tuple[list[str], tuple[tuple[str, str], ...]]:
    if not settings or not base_model or settings.worker_runner != "reinvent" or not smiles:
        return smiles, ()
    filename = "reinvent.prior" if base_model == "reinvent" else "mol2mol_scaffold_generic.prior"
    model_path = settings.reinvent_source_root / "priors" / filename
    return checkpoint_compatibility(settings, model_path, smiles)


def checkpoint_compatibility(
    settings: Settings, model_path: Path, smiles: list[str]
) -> tuple[list[str], tuple[tuple[str, str], ...]]:
    if settings.worker_runner != "reinvent" or not smiles:
        return smiles, ()
    checker = Path(__file__).resolve().parents[2] / "worker" / "smiles_compatibility.py"
    if not settings.reinvent_runtime_python.is_file() or not model_path.is_file():
        return smiles, ()
    environment = {
        "PYTHONPATH": os.pathsep.join((str(settings.reinvent_source_root), str(Path(__file__).resolve().parents[2]))),
        "PATH": os.environ.get("PATH", ""),
    }
    try:
        completed = subprocess.run(
            [str(settings.reinvent_runtime_python), str(checker), "--model", str(model_path)],
            input=json.dumps({"smiles": smiles}), capture_output=True, text=True,
            check=True, timeout=120, env=environment,
        )
        items = json.loads(completed.stdout)["items"]
    except (OSError, subprocess.SubprocessError, KeyError, json.JSONDecodeError) as error:
        raise AppError(
            status_code=422, code="model_compatibility_check_failed",
            message="REINVENTモデルとのSMILES互換性を確認できませんでした。",
        ) from error
    accepted = [str(item["smiles"]) for item in items if item.get("compatible")]
    rejected = tuple(
        (str(item["smiles"]), str(item.get("reason") or "model_incompatible"))
        for item in items if not item.get("compatible")
    )
    return accepted, rejected


def inspect_training_data(
    request: TrainingDataRequest, *, max_atoms: int, settings: Settings | None = None
) -> TrainingDataInspection:
    columns, values = _extract_values(request)
    if columns and not request.smiles_column:
        return TrainingDataInspection(columns, 0, [], 0, 0, [], standardization_counts([]))

    valid: list[str] = []
    seen: set[str] = set()
    results: list[StandardizationResult] = []
    duplicate_count = 0
    for value in values:
        result = standardize_smiles(value)
        results.append(result)
        if not result.accepted or not result.standardized_smiles:
            continue
        if result.standardized_smiles in seen:
            duplicate_count += 1
            continue
        seen.add(result.standardized_smiles)
        valid.append(result.standardized_smiles)
    valid, compatibility_rejections = _model_compatibility(settings, request.base_model, valid)
    compatibility_errors = tuple(reason for _, reason in compatibility_rejections)
    counts = standardization_counts(results)
    return TrainingDataInspection(
        columns=columns,
        total_rows=len(values),
        valid_smiles=valid,
        invalid_count=counts["invalid_count"] + counts["low_molecular_weight_excluded_count"],
        duplicate_count=duplicate_count,
        results=results,
        counts=counts,
        model_compatibility_error_count=len(compatibility_errors),
        compatibility_errors=compatibility_errors,
        compatibility_rejections=compatibility_rejections,
    )


def _next_model_name(settings: Settings, comment: str) -> str:
    date = datetime.now(ZoneInfo(settings.app_timezone)).strftime("%Y%m%d")
    base = f"{date}_{comment}"
    candidate = f"{base}.model"
    serial = 2
    existing = {path.name for path in settings.model_store_dir.glob(f"{base}*.model")}
    while candidate in existing:
        candidate = f"{base}_{serial}.model"
        serial += 1
    return candidate


def create_training_job(
    settings: Settings, request: TrainingJobCreateRequest
) -> JobRecord:
    inspection = inspect_training_data(request, max_atoms=settings.parent_max_atoms, settings=settings)
    if len(inspection.valid_smiles) < 2:
        raise AppError(
            status_code=422,
            code="insufficient_training_data",
            message="重複と不正構造を除いた有効なSMILESが2件以上必要です。",
        )

    shuffled = list(inspection.valid_smiles)
    random.Random(0).shuffle(shuffled)
    train_count = max(1, min(len(shuffled) - 1, round(len(shuffled) * request.train_fraction)))
    train_smiles = shuffled[:train_count]
    validation_smiles = shuffled[train_count:]

    settings.model_store_dir.mkdir(parents=True, exist_ok=True)
    model_name = _next_model_name(settings, request.comment)
    output_model = settings.model_store_dir / model_name
    output_model.touch(exist_ok=False)
    model_id = uuid4()
    job_id = uuid4()
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    for name in ("input", "config", "logs", "output"):
        (job_dir / name).mkdir(parents=True, exist_ok=False)
    (job_dir / "input" / "train.smi").write_text(
        "\n".join(train_smiles) + "\n", encoding="utf-8"
    )
    (job_dir / "input" / "validation.smi").write_text(
        "\n".join(validation_smiles) + "\n", encoding="utf-8"
    )
    accepted_path = job_dir / "input" / "standardization-accepted.csv"
    rejected_path = job_dir / "input" / "standardization-rejected.csv"
    compatible = set(inspection.valid_smiles)
    incompatibilities = dict(inspection.compatibility_rejections)
    audit_seen: set[str] = set()
    accepted_rows = []
    rejected_rows = []
    for row in inspection.results:
        if not row.accepted or not row.standardized_smiles:
            rejected_rows.append((row, row.reason or "invalid_smiles"))
        elif row.standardized_smiles in audit_seen:
            rejected_rows.append((row, "duplicate"))
        elif row.standardized_smiles not in compatible:
            rejected_rows.append((row, f"model_incompatible: {incompatibilities.get(row.standardized_smiles, 'unknown') }"))
            audit_seen.add(row.standardized_smiles)
        else:
            accepted_rows.append(row)
            audit_seen.add(row.standardized_smiles)
    with accepted_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Original_SMILES", "Standardized_SMILES", "MolecularWeight", "FormalCharge", "History"])
        writer.writerows(
            (row.original_smiles, row.standardized_smiles, row.molecular_weight, row.formal_charge, ";".join(row.history))
            for row in accepted_rows
        )
    with rejected_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Original_SMILES", "Standardized_SMILES", "MolecularWeight", "FormalCharge", "Reason", "History"])
        writer.writerows(
            (row.original_smiles, row.standardized_smiles or "", row.molecular_weight, row.formal_charge, reason, ";".join(row.history))
            for row, reason in rejected_rows
        )

    base_model_key = (
        "reinvent" if request.base_model == "reinvent" else "mol2mol_scaffold_generic"
    )
    conditions = {
        "base_model": request.base_model,
        "base_model_key": base_model_key,
        "model_name": model_name,
        "output_model_path": relative_to_repository(settings, output_model),
        "train_count": len(train_smiles),
        "validation_count": len(validation_smiles),
        "invalid_count": inspection.invalid_count,
        "duplicate_count": inspection.duplicate_count,
        "model_compatibility_error_count": inspection.model_compatibility_error_count,
        "standardization_counts": inspection.counts,
        "smiles_standardization": standardization_metadata(),
        "standardization_accepted_path": relative_to_repository(settings, accepted_path),
        "standardization_rejected_path": relative_to_repository(settings, rejected_path),
        "train_fraction": request.train_fraction,
        "num_epochs": request.num_epochs,
        "batch_size": request.batch_size,
        "learning_rate": request.learning_rate,
        "tanimoto_lower": request.tanimoto_lower,
        "tanimoto_upper": request.tanimoto_upper,
        "device": request.device,
        "smiles_column": request.smiles_column,
        "source_file_name": Path(request.file_name).name,
    }
    (job_dir / "input" / "request.json").write_text(
        json.dumps(
            {"comment": request.comment, "run_conditions": conditions},
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=request.comment,
        parent_smiles="",
        parent_smiles_canonical=None,
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, job_dir),
        job_type="transfer_learning",
        model_id=model_id,
    )
