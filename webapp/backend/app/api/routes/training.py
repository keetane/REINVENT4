from __future__ import annotations

from uuid import UUID

from fastapi import APIRouter, Request, status

from app.api.routes.jobs import _response
from app.db.jobs import JobRepository
from app.errors import AppError
from app.schemas.jobs import JobResponse
from app.schemas.training import (
    TrainingDataInspectionResponse,
    TrainingDataRequest,
    TrainingJobCreateRequest,
    TrainingProgressRequest,
)
from app.services.training import create_training_job, inspect_training_data
from app.services.smiles_standardization import standardization_metadata


router = APIRouter(tags=["transfer-learning"])


@router.post("/training-data/inspect", response_model=TrainingDataInspectionResponse)
def inspect_data(
    request: Request, payload: TrainingDataRequest
) -> TrainingDataInspectionResponse:
    result = inspect_training_data(
        payload, max_atoms=request.app.state.settings.parent_max_atoms,
        settings=request.app.state.settings,
    )
    return TrainingDataInspectionResponse(
        columns=result.columns,
        total_rows=result.total_rows,
        valid_count=len(result.valid_smiles),
        invalid_count=result.invalid_count,
        duplicate_count=result.duplicate_count,
        standardized_count=result.counts["standardized_count"],
        fragment_extracted_count=result.counts["fragment_extracted_count"],
        neutralized_count=result.counts["neutralized_count"],
        charged_structure_retained_count=result.counts["charged_structure_retained_count"],
        low_molecular_weight_excluded_count=result.counts["low_molecular_weight_excluded_count"],
        model_compatibility_error_count=result.model_compatibility_error_count,
        standardization_version=standardization_metadata()["version"],
    )


@router.post(
    "/training-jobs", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED
)
def post_training_job(
    request: Request, payload: TrainingJobCreateRequest
) -> JobResponse:
    return _response(create_training_job(request.app.state.settings, payload))


@router.post("/training-jobs/{job_id}/progress", response_model=JobResponse)
def report_training_progress(
    request: Request, job_id: UUID, payload: TrainingProgressRequest
) -> JobResponse:
    repository = JobRepository(request.app.state.settings)
    try:
        job = repository.get(job_id)
    except KeyError as error:
        raise AppError(
            status_code=404, code="job_not_found", message="学習ジョブが見つかりません。"
        ) from error
    if job.job_type != "transfer_learning":
        raise AppError(
            status_code=409,
            code="not_training_job",
            message="このジョブはTransfer Learningではありません。",
        )
    mean = payload.learning_mean
    return _response(
        repository.set_training_progress(
            job_id,
            epoch=payload.epoch,
            total_epochs=int(job.run_conditions["num_epochs"]),
            training_loss=mean.get("training"),
            validation_loss=mean.get("validation"),
        )
    )
