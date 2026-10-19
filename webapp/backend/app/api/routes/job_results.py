from __future__ import annotations

from collections.abc import Iterator
from datetime import datetime
import math
from pathlib import Path
from urllib.parse import quote
from uuid import UUID

from fastapi import APIRouter, Query, Request
from fastapi.responses import StreamingResponse

from app.db.jobs import JobRecord, JobRepository
from app.errors import AppError
from app.models.job import JobStatus
from app.schemas.results import (
    DescriptorSchema,
    MoleculeListResponse,
    MoleculeResponse,
    ResultSchemaResponse,
    StatisticsResponse,
)
from app.services.naming import confined_path, model_type_label, result_filename, slugify_comment
from app.services.result_query import (
    ResultDataset,
    ResultFilters,
    export_csv,
    export_sdf,
    load_result_dataset,
    molecule_svg,
    query_rows,
    scatter_for,
    statistics_for,
)
from app.services.results import ResultValidationError


router = APIRouter(prefix="/jobs/{job_id}/results", tags=["results"])
FILTER_QUERY_NAMES = {
    "MW": "mw",
    "RotB": "rotb",
    "HBD": "hbd",
    "HBA": "hba",
    "LogP": "logp",
    "TPSA": "tpsa",
    "QED": "qed",
    "SAScore": "sa_score",
    "FractionCSP3": "fraction_csp3",
    "RingCount": "ring_count",
    "AromaticRingCount": "aromatic_ring_count",
    "AliphaticRingCount": "aliphatic_ring_count",
    "HeavyAtomCount": "heavy_atom_count",
    "FormalCharge": "formal_charge",
    "LipinskiViolations": "lipinski_violations",
    "MolMR": "molmr",
    "ParentSimilarity": "parent_similarity",
    "ScaffoldSimilarity": "scaffold_similarity",
}


def _job_and_dataset(request: Request, job_id: UUID) -> tuple[JobRecord, ResultDataset]:
    repository = JobRepository(request.app.state.settings)
    try:
        job = repository.get(job_id)
    except KeyError as error:
        raise AppError(status_code=404, code="job_not_found", message="The requested job does not exist.") from error
    if job.status != JobStatus.COMPLETED or not job.result_csv_path:
        raise AppError(
            status_code=409,
            code="results_not_ready",
            message="Results are available only after the job completes.",
            context={"status": job.status.value},
        )
    try:
        path = confined_path(request.app.state.settings.repository_root, job.result_csv_path)
        parent_smiles = job.parent_smiles_canonical or job.parent_smiles_input
        return job, load_result_dataset(path, parent_smiles=parent_smiles)
    except (ValueError, ResultValidationError, UnicodeDecodeError) as error:
        raise AppError(
            status_code=500,
            code="results_unavailable",
            message="The completed result file could not be read.",
        ) from error


def _filters(
    *,
    request: Request,
    sort_by: str | None,
    sort_direction: str,
) -> ResultFilters:
    ranges = {}
    for descriptor, query_name in FILTER_QUERY_NAMES.items():
        ranges[descriptor] = (
            _query_number(request, f"{query_name}_min"),
            _query_number(request, f"{query_name}_max"),
        )
    return ResultFilters(
        ranges={name: bounds for name, bounds in ranges.items() if bounds != (None, None)},
        sort_by=sort_by,
        sort_direction=sort_direction,
    )


def _apply(dataset: ResultDataset, filters: ResultFilters):
    try:
        return query_rows(dataset, filters)
    except ValueError as error:
        raise AppError(status_code=422, code="invalid_result_query", message=str(error)) from error


def _query_number(request: Request, name: str) -> float | None:
    raw = request.query_params.get(name)
    if raw in {None, ""}:
        return None
    try:
        value = float(raw)
    except ValueError as error:
        raise AppError(
            status_code=422,
            code="invalid_result_query",
            message=f"{name} must be a number.",
        ) from error
    if not math.isfinite(value):
        raise AppError(
            status_code=422,
            code="invalid_result_query",
            message="Filter bounds must be finite numbers.",
        )
    return value


def _route_filters(values: dict) -> ResultFilters:
    return _filters(
        request=values["request"],
        sort_by=values["sort_by"],
        sort_direction=values["sort_direction"],
    )


@router.get("/schema", response_model=ResultSchemaResponse)
def get_schema(request: Request, job_id: UUID) -> ResultSchemaResponse:
    _, dataset = _job_and_dataset(request, job_id)
    return ResultSchemaResponse(
        job_id=job_id,
        columns=list(dataset.columns),
        descriptors=[
            DescriptorSchema(
                name=item.name,
                source=item.source,
                aliases=list(item.aliases),
                minimum=item.minimum,
                maximum=item.maximum,
                missing_count=item.missing_count,
            )
            for item in dataset.descriptors
        ],
        total_rows=dataset.total_rows,
        valid_rows=len(dataset.rows),
        invalid_smiles=dataset.invalid_smiles,
        warnings=list(dataset.warnings),
        rdkit_version=__import__("rdkit").__version__,
    )


@router.get("/molecules", response_model=MoleculeListResponse)
def get_molecules(
    request: Request,
    job_id: UUID,
    limit: int = Query(default=50, ge=1, le=500),
    offset: int = Query(default=0, ge=0),
    sort_by: str | None = None,
    sort_direction: str = "asc",
    mw_min: float | None = None,
    mw_max: float | None = None,
    rotb_min: float | None = None,
    rotb_max: float | None = None,
    hbd_min: float | None = None,
    hbd_max: float | None = None,
    hba_min: float | None = None,
    hba_max: float | None = None,
    logp_min: float | None = None,
    logp_max: float | None = None,
    tpsa_min: float | None = None,
    tpsa_max: float | None = None,
) -> MoleculeListResponse:
    _, dataset = _job_and_dataset(request, job_id)
    filters = _route_filters(locals())
    rows = _apply(dataset, filters)
    return MoleculeListResponse(
        job_id=job_id,
        items=[
            MoleculeResponse(
                row_number=row.row_number,
                smiles=row.smiles,
                canonical_smiles=row.canonical_smiles,
                structure_svg=molecule_svg(row.molecule),
                values=row.values,
            )
            for row in rows[offset : offset + limit]
        ],
        total=len(rows),
        limit=limit,
        offset=offset,
        invalid_smiles=dataset.invalid_smiles,
    )


@router.get("/statistics", response_model=StatisticsResponse)
def get_statistics(
    request: Request,
    job_id: UUID,
    sort_by: str | None = None,
    sort_direction: str = "asc",
    mw_min: float | None = None,
    mw_max: float | None = None,
    rotb_min: float | None = None,
    rotb_max: float | None = None,
    hbd_min: float | None = None,
    hbd_max: float | None = None,
    hba_min: float | None = None,
    hba_max: float | None = None,
    logp_min: float | None = None,
    logp_max: float | None = None,
    tpsa_min: float | None = None,
    tpsa_max: float | None = None,
    scatter_x: str = "MW",
    scatter_y: str = "LogP",
) -> StatisticsResponse:
    _, dataset = _job_and_dataset(request, job_id)
    filters = _route_filters(locals())
    rows = _apply(dataset, filters)
    summaries, histograms = statistics_for(dataset, rows)
    try:
        scatter, scatter_sampled = scatter_for(dataset, rows, scatter_x, scatter_y)
    except ValueError as error:
        raise AppError(status_code=422, code="invalid_result_query", message=str(error)) from error
    return StatisticsResponse(
        job_id=job_id,
        total_rows=dataset.total_rows,
        filtered_rows=len(rows),
        invalid_smiles=dataset.invalid_smiles,
        summaries=summaries,
        histograms=histograms,
        scatter_x=scatter_x,
        scatter_y=scatter_y,
        scatter_sampled=scatter_sampled,
        scatter=scatter,
    )


def _chunks(value: str, size: int = 64 * 1024) -> Iterator[str]:
    for start in range(0, len(value), size):
        yield value[start : start + size]


def _download_filename(request: Request, job: JobRecord) -> str:
    """Build the current result name even for jobs created before the naming update."""
    try:
        created_at = datetime.fromisoformat(job.created_at)
        model_key = str(job.run_conditions.get("model_key", "reinvent"))
        return result_filename(
            comment_slug=job.comment_slug,
            created_at=created_at,
            timezone_name=request.app.state.settings.app_timezone,
            tab_type="MolGen",
            model_type=model_type_label(model_key),
        )
    except (TypeError, ValueError, KeyError):
        return job.result_csv_path.rsplit("/", 1)[-1]


def _content_disposition(filename: str) -> str:
    """Return a browser-compatible attachment name with a UTF-8 fallback."""
    path = Path(filename)
    ascii_stem = slugify_comment(path.stem, ascii_only=True, max_length=180)
    ascii_filename = f"{ascii_stem}{path.suffix}"
    return (
        f'attachment; filename="{ascii_filename}"; '
        f"filename*=UTF-8''{quote(filename)}"
    )


@router.get("/download.csv")
def download_csv(
    request: Request,
    job_id: UUID,
    sort_by: str | None = None,
    sort_direction: str = "asc",
    mw_min: float | None = None,
    mw_max: float | None = None,
    rotb_min: float | None = None,
    rotb_max: float | None = None,
    hbd_min: float | None = None,
    hbd_max: float | None = None,
    hba_min: float | None = None,
    hba_max: float | None = None,
    logp_min: float | None = None,
    logp_max: float | None = None,
    tpsa_min: float | None = None,
    tpsa_max: float | None = None,
):
    job, dataset = _job_and_dataset(request, job_id)
    filters = _route_filters(locals())
    rows = _apply(dataset, filters)
    content = export_csv(dataset.columns, rows, (item.name for item in dataset.descriptors))
    filename = _download_filename(request, job)
    return StreamingResponse(
        _chunks(content),
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": _content_disposition(filename)},
    )


@router.get("/download.sdf")
def download_sdf(
    request: Request,
    job_id: UUID,
    sort_by: str | None = None,
    sort_direction: str = "asc",
    mw_min: float | None = None,
    mw_max: float | None = None,
    rotb_min: float | None = None,
    rotb_max: float | None = None,
    hbd_min: float | None = None,
    hbd_max: float | None = None,
    hba_min: float | None = None,
    hba_max: float | None = None,
    logp_min: float | None = None,
    logp_max: float | None = None,
    tpsa_min: float | None = None,
    tpsa_max: float | None = None,
):
    job, dataset = _job_and_dataset(request, job_id)
    filters = _route_filters(locals())
    rows = _apply(dataset, filters)
    content = export_sdf(rows)
    filename = _download_filename(request, job).removesuffix(".csv") + ".sdf"
    return StreamingResponse(
        _chunks(content),
        media_type="chemical/x-mdl-sdfile",
        headers={"Content-Disposition": _content_disposition(filename)},
    )
