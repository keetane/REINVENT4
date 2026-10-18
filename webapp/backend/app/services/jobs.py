from __future__ import annotations

from uuid import UUID, uuid4

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.db.models import ModelRepository
from app.errors import AppError
from app.reinvent_adapter.models import MODEL_SPECS
from app.schemas.jobs import JobCreateRequest
from app.services.naming import slugify_comment
from app.services.parents import validate_parent_smiles
from app.services.workspaces import create_job_workspace, relative_to_repository


def create_job(settings: Settings, request: JobCreateRequest) -> JobRecord:
    conditions = dict(request.run_conditions)
    requested_model_id = str(
        conditions.get("model_id") or f"standard:{conditions.get('model_key', 'mol2mol_high')}"
    )
    conditions.pop("resolved_model_path", None)
    custom_model_id: UUID | None = None
    if requested_model_id.startswith("standard:"):
        model_key = requested_model_id.removeprefix("standard:")
        try:
            spec = MODEL_SPECS[model_key]
        except KeyError as error:
            raise AppError(
                status_code=422,
                code="unknown_model",
                message="選択した標準モデルは利用できません。",
            ) from error
        model_type = spec.mode
        resolved_path = spec.relative_path
    else:
        try:
            custom_model_id = UUID(requested_model_id)
            model = ModelRepository(settings).get(custom_model_id)
        except (ValueError, KeyError) as error:
            raise AppError(
                status_code=404,
                code="model_not_found",
                message="選択した学習済みモデルが見つかりません。",
            ) from error
        model_key = model.base_model_key
        model_type = model.model_type
        resolved_path = model.file_path

    conditions.update(
        {
            "model_id": requested_model_id,
            "model_key": model_key,
            "model_type": model_type,
            "resolved_model_path": resolved_path,
        }
    )
    if model_type == "reinvent":
        parent_input = ""
        parent_canonical = None
    else:
        validation_key = model_key if model_key in {"libinvent", "linkinvent"} else "mol2mol_high"
        parent = validate_parent_smiles(
            request.parent_smiles,
            max_atoms=settings.parent_max_atoms,
            model_key=validation_key,
        )
        parent_input = parent.input_smiles
        parent_canonical = parent.canonical_smiles
    job_id = uuid4()
    slug = slugify_comment(request.comment)
    workspace = create_job_workspace(
        settings,
        job_id=job_id,
        parent_smiles=parent_canonical or "",
        comment=request.comment,
        run_conditions=conditions,
    )
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=slug,
        parent_smiles=parent_input,
        parent_smiles_canonical=parent_canonical,
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, workspace),
        model_id=custom_model_id,
    )
