from __future__ import annotations

from uuid import UUID

from fastapi import APIRouter, Request

from app.db.models import ModelRepository
from app.errors import AppError
from app.reinvent_adapter.models import MODEL_SPECS
from app.schemas.models import ModelDeleteResponse, ModelListResponse, ModelResponse
from app.services.naming import confined_path


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
