from __future__ import annotations

from datetime import datetime
from uuid import UUID, uuid4
from zoneinfo import ZoneInfo

from fastapi import APIRouter, Query, Request
from fastapi.responses import FileResponse

from app.db.jobs import JobRepository
from app.db.models import ModelRepository
from app.errors import AppError
from app.reinvent_adapter.models import MODEL_SPECS
from app.schemas.models import ModelDeleteResponse, ModelListResponse, ModelResponse
from app.services.naming import confined_path, infer_model_metadata, next_artifact_filename, slugify_comment
from app.services.workspaces import relative_to_repository


router = APIRouter(prefix="/models", tags=["models"])

MODEL_LABELS = {
    "mol2mol_high": ("Mol2Mol high similarity", "Parent分子に近い構造を生成します。"),
    "mol2mol_medium": ("Mol2Mol medium", "Parent分子から中程度の変換を行います。"),
    "mol2mol_mmp": ("Mol2Mol MMP", "Matched molecular pair型の変換を行います。"),
    "mol2mol_scaffold": ("Mol2Mol scaffold", "Parent分子の骨格変換を行います。"),
    "mol2mol_scaffold_generic": ("Mol2Mol scaffold generic", "Transfer LearningとRLに対応する汎用Mol2Molモデルです。"),
    "reinvent": ("REINVENT Prior", "Parent不要でPriorから分子を生成します。"),
    "libinvent": ("LibInvent", "骨格を固定して側鎖を生成します。"),
    "linkinvent": ("LinkInvent", "側鎖を固定してリンカー骨格を生成します。"),
}
MAX_MODEL_UPLOAD_BYTES = 2 * 1024 * 1024 * 1024


@router.get("", response_model=ModelListResponse)
def list_models(request: Request) -> ModelListResponse:
    standard = []
    for key in (
        "mol2mol_high",
        "mol2mol_medium",
        "mol2mol_mmp",
        "mol2mol_scaffold",
        "mol2mol_scaffold_generic",
        "reinvent",
        "libinvent",
        "linkinvent",
    ):
        label, description = MODEL_LABELS[key]
        standard.append(
            ModelResponse(
                id=f"standard:{key}",
                name=label,
                model_type=MODEL_SPECS[key].mode,
                source="standard",
                description=description,
                base_model_key=key,
            )
        )
    custom = [
        ModelResponse(
            id=str(model.id),
            name=model.name,
            model_type=model.model_type,
            source=("reinforcement_learning" if model.training_conditions.get("source") == "reinforcement_learning" else "transfer_learning"),
            description=(
                f"{model.model_type.upper()} Reinforcement Learningモデル"
                if model.training_conditions.get("source") == "reinforcement_learning"
                else f"{model.model_type.upper()} Transfer Learningモデル"
            ),
            created_at=model.created_at,
            training_data_count=model.training_data_count,
            validation_data_count=model.validation_data_count,
            training_conditions=model.training_conditions,
            base_model_key=model.base_model_key,
            deletable=True,
        )
        for model in ModelRepository(request.app.state.settings).list()
    ]
    return ModelListResponse(items=standard + custom)


def _safe_upload_comment(value: str) -> str:
    cleaned = slugify_comment(value)
    if cleaned == "result" and not value.strip():
        raise AppError(status_code=422, code="invalid_model_name", message="モデル名を入力してください。")
    return cleaned


@router.post("/upload", response_model=ModelResponse, status_code=201)
async def upload_model(
    request: Request,
    name: str | None = Query(default=None, min_length=1, max_length=100),
    model_type: str | None = Query(default=None, pattern="^(reinvent|mol2mol)$"),
    base_model_key: str | None = Query(default=None),
    source: str | None = Query(default="transfer_learning", pattern="^(transfer_learning|reinforcement_learning)$"),
) -> ModelResponse:
    """Import a compatible REINVENT checkpoint into the model registry.

    The body is the raw model file so the endpoint does not require the
    optional multipart package and can handle large checkpoint files directly.
    """
    settings = request.app.state.settings
    original_filename = request.headers.get("x-model-filename", "uploaded-model")
    inferred = infer_model_metadata(original_filename)
    if inferred is not None:
        name = inferred["comment"]
        model_type = inferred["model_type"]
        base_model_key = inferred["base_model_key"]
        source = inferred["source"]
    if not name or not model_type or not base_model_key or not source:
        raise AppError(
            status_code=422,
            code="model_metadata_required",
            message="ファイル名から用途・種別を判定できないため、モデル情報を指定してください。",
        )
    spec = MODEL_SPECS.get(base_model_key)
    if spec is None or spec.mode != model_type:
        raise AppError(status_code=422, code="incompatible_model_metadata", message="モデル種別とBase Modelの組み合わせが不正です。")
    content_length = request.headers.get("content-length")
    if content_length:
        try:
            too_large = int(content_length) > MAX_MODEL_UPLOAD_BYTES
        except ValueError as error:
            raise AppError(status_code=400, code="invalid_content_length", message="アップロードサイズを確認できません。") from error
        if too_large:
            raise AppError(status_code=413, code="model_upload_too_large", message="モデルファイルは2 GiB以下にしてください。")
    data = await request.body()
    if not data:
        raise AppError(status_code=422, code="empty_model_upload", message="モデルファイルが空です。")
    if len(data) > MAX_MODEL_UPLOAD_BYTES:
        raise AppError(status_code=413, code="model_upload_too_large", message="モデルファイルは2 GiB以下にしてください。")

    comment = _safe_upload_comment(name)
    extension = ".chkpt" if source == "reinforcement_learning" else ".model"
    settings.model_store_dir.mkdir(parents=True, exist_ok=True)
    if comment in {item.name for item in ModelRepository(settings).list()}:
        raise AppError(status_code=409, code="duplicate_model_name", message="同じ表示名のモデルが既に登録されています。")
    candidate = next_artifact_filename(
        settings.model_store_dir,
        tab_type="RL" if source == "reinforcement_learning" else "TL",
        model_key=base_model_key,
        comment=comment,
        date=datetime.now(ZoneInfo(settings.app_timezone)),
        extension=extension,
    )
    destination = settings.model_store_dir / candidate
    destination.write_bytes(data)

    model_id = uuid4()
    job_id = uuid4()
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    job_dir.mkdir(parents=True, exist_ok=False)
    conditions = {
        "source": source,
        "model_upload": True,
        "model_name": comment,
        "model_file_name": candidate,
        "base_model": model_type,
        "base_model_key": base_model_key,
        "original_filename": original_filename,
    }
    JobRepository(settings).create_completed_import(
        job_id=job_id,
        comment=comment,
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, job_dir),
        model_id=model_id,
    )
    model = ModelRepository(settings).register(
        model_id=model_id,
        name=comment,
        model_type=model_type,
        training_data_count=1,
        validation_data_count=0,
        training_conditions=conditions,
        file_path=relative_to_repository(settings, destination),
        config_path="",
        log_path="",
        source_job_id=job_id,
        base_model_key=base_model_key,
    )
    return ModelResponse(
        id=str(model.id), name=model.name, model_type=model.model_type,
        source=source, description=f"アップロードされた{model_type.upper()}モデルです。",
        created_at=model.created_at, training_data_count=model.training_data_count,
        validation_data_count=model.validation_data_count, training_conditions=model.training_conditions,
        base_model_key=model.base_model_key, deletable=True,
    )


@router.get("/{model_id}/artifacts/{kind}")
def download_model_artifact(request: Request, model_id: UUID, kind: str) -> FileResponse:
    if kind != "model":
        raise AppError(status_code=404, code="model_artifact_not_found", message="指定した成果物はありません。")
    repository = ModelRepository(request.app.state.settings)
    try:
        model = repository.get(model_id)
    except KeyError as error:
        raise AppError(status_code=404, code="model_not_found", message="モデルが見つかりません。") from error
    path = confined_path(request.app.state.settings.repository_root, model.file_path)
    if not path.is_file() or path.is_symlink() or not path.is_relative_to(request.app.state.settings.model_store_dir):
        raise AppError(status_code=404, code="model_artifact_not_found", message="モデルファイルが見つかりません。")
    return FileResponse(path, filename=path.name, media_type="application/octet-stream")


@router.delete("/{model_id}", response_model=ModelDeleteResponse)
def delete_model(request: Request, model_id: UUID) -> ModelDeleteResponse:
    settings = request.app.state.settings
    repository = ModelRepository(settings)
    try:
        model = repository.get(model_id)
    except KeyError as error:
        raise AppError(
            status_code=404, code="model_not_found", message="モデルが見つかりません。"
        ) from error
    if repository.is_in_use(model_id):
        raise AppError(
            status_code=409,
            code="model_in_use",
            message="学習中または生成ジョブで使用中のモデルは削除できません。",
        )
    model_path = confined_path(settings.repository_root, model.file_path)
    if not model_path.is_relative_to(settings.model_store_dir):
        raise AppError(
            status_code=500,
            code="invalid_model_path",
            message="モデルファイルの保存先が不正です。",
        )
    if model_path.is_file() and not model_path.is_symlink():
        model_path.unlink()
    repository.delete(model_id)
    return ModelDeleteResponse(id=model_id, deleted=True)
