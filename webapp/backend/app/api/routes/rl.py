from __future__ import annotations

from uuid import UUID

from fastapi import APIRouter, Request, status

from app.api.routes.jobs import _response
from app.db.jobs import JobRepository
from app.errors import AppError
from app.schemas.jobs import JobResponse
from app.schemas.rl import (
    QsarPredictionRequest,
    QsarPredictionResponse,
    RlJobCreateRequest,
    RlPreviewRequest,
    RlPreviewResponse,
    RlProgressRequest,
)
from app.services.rl import available_torch_devices, create_rl_job, predict_qsar, preview_rl_scoring


router = APIRouter(tags=["reinforcement-learning"])


@router.get("/rl/capabilities")
def rl_capabilities(request: Request) -> dict:
    devices = available_torch_devices(request.app.state.settings)
    return {
        "devices": [
            {"key": key, "label": key.upper(), "available": available}
            for key, available in devices.items()
        ]
    }


@router.post("/rl-jobs", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED)
def post_rl_job(request: Request, payload: RlJobCreateRequest) -> JobResponse:
    return _response(create_rl_job(request.app.state.settings, payload))


@router.post("/rl-jobs/{job_id}/progress", response_model=JobResponse)
def report_rl_progress(
    request: Request, job_id: UUID, payload: RlProgressRequest
) -> JobResponse:
    repository = JobRepository(request.app.state.settings)
    try:
        job = repository.get(job_id)
    except KeyError as error:
        raise AppError(status_code=404, code="job_not_found", message="RLジョブが見つかりません。") from error
    if job.job_type != "reinforcement_learning":
        raise AppError(status_code=409, code="not_rl_job", message="このジョブはReinforcement Learningではありません。")
    metrics = {
        "stage": payload.stage,
        "components": payload.components,
        "learning": payload.learning,
        "fraction_valid_smiles": payload.fraction_valid_smiles,
        "collected_smiles_in_memory": payload.collected_smiles_in_memory,
        "mean_total_score": payload.components.get("total_score"),
    }
    return _response(
        repository.set_rl_progress(
            job_id,
            step=payload.step,
            maximum_steps=int(job.run_conditions["maximum_steps"]),
            metrics=metrics,
        )
    )


@router.post("/qsar/predict", response_model=QsarPredictionResponse)
def qsar_predict(request: Request, payload: QsarPredictionRequest) -> QsarPredictionResponse:
    return QsarPredictionResponse(
        model_id=payload.model_id,
        predictions=predict_qsar(request.app.state.settings, payload),
    )


@router.post("/rl/preview", response_model=RlPreviewResponse)
def rl_preview(request: Request, payload: RlPreviewRequest) -> RlPreviewResponse:
    return RlPreviewResponse(rows=preview_rl_scoring(request.app.state.settings, payload))
