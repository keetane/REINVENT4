from __future__ import annotations

from uuid import uuid4

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.schemas.jobs import JobCreateRequest
from app.services.naming import slugify_comment
from app.services.parents import validate_parent_smiles
from app.services.workspaces import create_job_workspace, relative_to_repository


def create_job(settings: Settings, request: JobCreateRequest) -> JobRecord:
    model_key = str(request.run_conditions.get("model_key", "mol2mol_high"))
    parent = validate_parent_smiles(
        request.parent_smiles,
        max_atoms=settings.parent_max_atoms,
        model_key=model_key,
    )
    job_id = uuid4()
    slug = slugify_comment(request.comment)
    workspace = create_job_workspace(
        settings,
        job_id=job_id,
        parent_smiles=parent.canonical_smiles,
        comment=request.comment,
        run_conditions=request.run_conditions,
    )
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=slug,
        parent_smiles=parent.input_smiles,
        parent_smiles_canonical=parent.canonical_smiles,
        run_conditions=request.run_conditions,
        job_dir=relative_to_repository(settings, workspace),
    )
