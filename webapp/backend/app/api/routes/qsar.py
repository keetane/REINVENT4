from __future__ import annotations

import json
import shutil
import subprocess
from functools import lru_cache
from uuid import UUID, uuid4

from fastapi import APIRouter, Request, status
from fastapi.responses import FileResponse

from app.api.routes.jobs import _response
from app.db.jobs import JobRepository
from app.db.qsar_models import QsarModelRecord, QsarModelRepository
from app.errors import AppError
from app.schemas.jobs import JobResponse
from app.schemas.qsar import (
    QsarCapabilitiesResponse,
    QsarCapability,
    QsarComparisonCreateRequest,
    QsarComparisonRegistrationRequest,
    QsarDataInspectionResponse,
    QsarDataRequest,
    QsarJobCreateRequest,
    QsarModelDeleteResponse,
    QsarModelListResponse,
    QsarModelResponse,
)
from app.services.naming import confined_path
from app.services.qsar import _next_model_name, create_qsar_comparison_job, create_qsar_job, process_qsar_data
from app.services.smiles_standardization import standardization_metadata
from app.services.workspaces import relative_to_repository


router = APIRouter(tags=["qsar"])


@lru_cache(maxsize=4)
def _runtime_versions(python: str) -> dict[str, str | None]:
    # PackageNotFoundError would abort the compact command, so query independently.
    result: dict[str, str | None] = {}
    for package in ("scikit-learn", "lightgbm", "xgboost", "catboost", "chemprop", "qsartuna"):
        query = (
            "import importlib.metadata as m;"
            f"print(m.version({package!r}))"
        )
        try:
            completed = subprocess.run(
                [python, "-c", query],
                check=True,
                capture_output=True,
                text=True,
                timeout=10,
            )
            result[package] = completed.stdout.strip()
        except (OSError, subprocess.SubprocessError):
            result[package] = None
    return result


def _model_response(repository: QsarModelRepository, model: QsarModelRecord) -> QsarModelResponse:
    return QsarModelResponse(
        id=model.id,
        name=model.name,
        algorithm=model.algorithm,
        task_type=model.task_type,
        target_name=model.target_name,
        target_unit=model.target_unit,
        created_at=model.created_at,
        training_data_count=model.training_data_count,
        test_data_count=model.test_data_count,
        metrics=model.metrics,
        metadata=model.metadata,
        ad_threshold=model.ad_threshold,
        in_use=repository.is_in_use(model.id),
    )


@router.get("/qsar/capabilities", response_model=QsarCapabilitiesResponse)
def capabilities(request: Request) -> QsarCapabilitiesResponse:
    settings = request.app.state.settings
    versions = _runtime_versions(str(settings.reinvent_runtime_python))
    chemprop_versions = _runtime_versions(str(settings.chemprop_runtime_python)) if settings.chemprop_runtime_python and settings.chemprop_runtime_python.is_file() else {}
    qsartuna_versions = _runtime_versions(str(settings.qsartuna_runtime_python)) if settings.qsartuna_runtime_python and settings.qsartuna_runtime_python.is_file() else {}
    sklearn_version = versions.get("scikit-learn")
    chemprop_version = chemprop_versions.get("chemprop")
    chemprop_available = bool(chemprop_version and int(chemprop_version.split(".", 1)[0]) >= 2)
    qsartuna_version = qsartuna_versions.get("qsartuna")
    algorithms = [
        QsarCapability(key="random_forest", label="Random Forest", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(key="svm", label="SVM", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(key="gradient_boosting", label="Gradient Boosting", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(key="extra_trees", label="Extra Trees", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(key="lightgbm", label="LightGBM", available=bool(versions.get("lightgbm")), version=versions.get("lightgbm"), reason=None if versions.get("lightgbm") else "LightGBMがインストールされていません。"),
        QsarCapability(key="xgboost", label="XGBoost", available=bool(versions.get("xgboost")), version=versions.get("xgboost"), reason=None if versions.get("xgboost") else "XGBoostがインストールされていません。"),
        QsarCapability(key="catboost", label="CatBoost", available=bool(versions.get("catboost")), version=versions.get("catboost"), reason=None if versions.get("catboost") else "CatBoostがインストールされていません。"),
        QsarCapability(key="kernel_ridge", label="Kernel Ridge", available=bool(sklearn_version), version=sklearn_version, task_types=["regression"]),
        QsarCapability(key="gaussian_process", label="Gaussian Process", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(key="ridge", label="Ridge", available=bool(sklearn_version), version=sklearn_version, task_types=["regression"]),
        QsarCapability(key="elastic_net", label="Elastic Net", available=bool(sklearn_version), version=sklearn_version, task_types=["regression"]),
        QsarCapability(key="logistic_regression", label="Logistic Regression", available=bool(sklearn_version), version=sklearn_version, task_types=["classification"]),
        QsarCapability(key="mlp", label="MLP", available=bool(sklearn_version), version=sklearn_version),
        QsarCapability(
            key="chemprop_v2",
            label="Chemprop v2",
            available=chemprop_available,
            version=chemprop_version,
            reason=(
                None if chemprop_available else f"検出したChemprop {chemprop_version}はv2ではありません。"
                if chemprop_version else "Chemprop v2の独立環境がインストールされていません。"
            ),
        ),
        QsarCapability(
            key="qsartuna",
            label="QSARtuna",
            available=bool(qsartuna_version),
            version=qsartuna_version,
            reason=None if qsartuna_version else "QSARtunaの独立環境がインストールされていません。",
        ),
    ]
    return QsarCapabilitiesResponse(
        algorithms=algorithms,
        devices=[QsarCapability(key="cpu", label="CPU", available=True, reason="scikit-learnモデルはCPUで学習します。")],
    )


@router.post("/qsar-data/inspect", response_model=QsarDataInspectionResponse)
def inspect_qsar_data(request: Request, payload: QsarDataRequest) -> QsarDataInspectionResponse:
    result = process_qsar_data(payload, max_atoms=request.app.state.settings.parent_max_atoms)
    return QsarDataInspectionResponse(
        columns=result.columns,
        total_rows=result.total_rows,
        valid_count=len(result.rows),
        invalid_smiles_count=result.invalid_smiles_count,
        missing_target_count=result.missing_target_count,
        invalid_target_count=result.invalid_target_count,
        duplicate_count=result.duplicate_count,
        conflicting_duplicate_count=result.conflicting_duplicate_count,
        standardized_count=result.standardization_counts["standardized_count"],
        fragment_extracted_count=result.standardization_counts["fragment_extracted_count"],
        neutralized_count=result.standardization_counts["neutralized_count"],
        charged_structure_retained_count=result.standardization_counts["charged_structure_retained_count"],
        low_molecular_weight_excluded_count=result.standardization_counts["low_molecular_weight_excluded_count"],
        final_accepted_count=len(result.rows),
        standardization_version=standardization_metadata()["version"],
        preview=result.preview,
        rejected_rows=result.rejected_rows,
    )


@router.post("/qsar-jobs", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED)
def post_qsar_job(request: Request, payload: QsarJobCreateRequest) -> JobResponse:
    return _response(create_qsar_job(request.app.state.settings, payload))


@router.post("/qsar-comparison-jobs", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED)
def post_qsar_comparison_job(request: Request, payload: QsarComparisonCreateRequest) -> JobResponse:
    return _response(create_qsar_comparison_job(request.app.state.settings, payload))


def _comparison(request: Request, job_id: UUID):
    try:
        job = JobRepository(request.app.state.settings).get(job_id)
    except KeyError as error:
        raise AppError(status_code=404, code="qsar_comparison_not_found", message="モデル比較ジョブが見つかりません。") from error
    if job.job_type != "qsar_comparison":
        raise AppError(status_code=404, code="qsar_comparison_not_found", message="モデル比較ジョブが見つかりません。")
    if job.status.value != "completed":
        raise AppError(status_code=409, code="qsar_comparison_not_completed", message="モデル比較はまだ完了していません。")
    summary_path = confined_path(request.app.state.settings.repository_root, job.run_conditions["summary_path"])
    if not summary_path.is_file():
        raise AppError(status_code=404, code="qsar_comparison_artifact_not_found", message="比較結果が見つかりません。")
    return job, json.loads(summary_path.read_text(encoding="utf-8"))


@router.get("/qsar-comparisons/{job_id}")
def get_qsar_comparison(request: Request, job_id: UUID) -> dict:
    _, summary = _comparison(request, job_id)
    return summary


@router.post("/qsar-comparisons/{job_id}/register", response_model=QsarModelResponse)
def register_qsar_comparison_model(request: Request, job_id: UUID, payload: QsarComparisonRegistrationRequest) -> QsarModelResponse:
    job, summary = _comparison(request, job_id)
    result = next((item for item in summary["results"] if item["algorithm"] == payload.algorithm and item["status"] == "completed"), None)
    if result is None:
        raise AppError(status_code=422, code="qsar_comparison_model_unavailable", message="選択したモデルは正常に学習されていません。")
    settings = request.app.state.settings
    source_model = confined_path(settings.repository_root, result["model_path"])
    evaluation_path = confined_path(settings.repository_root, result["evaluation_path"])
    predictions_path = confined_path(settings.repository_root, result["predictions_path"])
    if not source_model.is_file() or source_model.stat().st_size == 0:
        raise AppError(status_code=404, code="qsar_comparison_artifact_not_found", message="選択したモデルファイルが見つかりません。")
    settings.qsar_model_store_dir.mkdir(parents=True, exist_ok=True)
    name = _next_model_name(settings, payload.comment)
    destination = settings.qsar_model_store_dir / f"{name}.joblib"
    shutil.copy2(source_model, destination)
    evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
    conditions = job.run_conditions
    model_id = uuid4()
    repository = QsarModelRepository(settings)
    model = repository.register(
        model_id=model_id,
        name=name,
        algorithm=str(payload.algorithm),
        task_type=str(conditions["task_type"]),
        target_name=str(conditions["target_name"]),
        target_unit=str(conditions["target_unit"]),
        training_data_count=int(evaluation["metrics"]["training_data_count"]),
        test_data_count=int(evaluation["metrics"]["test_data_count"]),
        metrics=evaluation["metrics"],
        metadata={**evaluation["metadata"], "comparison_job_id": str(job_id), "model_id": str(model_id)},
        file_path=relative_to_repository(settings, destination),
        dataset_path=str(conditions["dataset_path"]),
        predictions_path=relative_to_repository(settings, predictions_path),
        evaluation_path=relative_to_repository(settings, evaluation_path),
        source_job_id=job_id,
        ad_threshold=float(conditions["ad_threshold"]),
    )
    return _model_response(repository, model)


@router.get("/qsar-models", response_model=QsarModelListResponse)
def list_qsar_models(request: Request) -> QsarModelListResponse:
    repository = QsarModelRepository(request.app.state.settings)
    return QsarModelListResponse(items=[_model_response(repository, model) for model in repository.list()])


def _get_model(request: Request, model_id: UUID) -> tuple[QsarModelRepository, QsarModelRecord]:
    repository = QsarModelRepository(request.app.state.settings)
    try:
        return repository, repository.get(model_id)
    except KeyError as error:
        raise AppError(status_code=404, code="qsar_model_not_found", message="QSARモデルが見つかりません。") from error


@router.get("/qsar-models/{model_id}", response_model=QsarModelResponse)
def get_qsar_model(request: Request, model_id: UUID) -> QsarModelResponse:
    repository, model = _get_model(request, model_id)
    return _model_response(repository, model)


def _artifact(request: Request, model_id: UUID, kind: str) -> FileResponse:
    _, model = _get_model(request, model_id)
    attribute = {
        "model": "file_path",
        "dataset": "dataset_path",
        "predictions": "predictions_path",
        "evaluation": "evaluation_path",
    }.get(kind)
    if attribute is None:
        raise AppError(status_code=404, code="qsar_artifact_not_found", message="指定した成果物はありません。")
    path = confined_path(request.app.state.settings.repository_root, getattr(model, attribute))
    if not path.is_file() or path.is_symlink():
        raise AppError(status_code=404, code="qsar_artifact_not_found", message="成果物ファイルが見つかりません。")
    media_type = "application/octet-stream" if kind == "model" else "text/csv" if path.suffix == ".csv" else "application/json"
    return FileResponse(path, filename=path.name, media_type=media_type)


@router.get("/qsar-models/{model_id}/artifacts/{kind}")
def download_qsar_artifact(request: Request, model_id: UUID, kind: str) -> FileResponse:
    return _artifact(request, model_id, kind)


@router.delete("/qsar-models/{model_id}", response_model=QsarModelDeleteResponse)
def delete_qsar_model(request: Request, model_id: UUID) -> QsarModelDeleteResponse:
    repository, model = _get_model(request, model_id)
    if repository.is_in_use(model_id):
        raise AppError(status_code=409, code="qsar_model_in_use", message="強化学習ジョブで使用中のQSARモデルは削除できません。")
    settings = request.app.state.settings
    model_path = confined_path(settings.repository_root, model.file_path)
    if not model_path.is_relative_to(settings.qsar_model_store_dir):
        raise AppError(status_code=500, code="invalid_qsar_model_path", message="QSARモデルの保存先が不正です。")
    if model_path.is_file() and not model_path.is_symlink():
        model_path.unlink()
    repository.delete(model_id)
    return QsarModelDeleteResponse(id=model_id, deleted=True)
