from pathlib import Path

from fastapi.testclient import TestClient

from app.db.jobs import JobRepository
from app.main import create_app
from app.models.job import JobStatus
from app.schemas.jobs import JobCreateRequest
from app.services.jobs import create_job
from app.services.workspaces import relative_to_repository


def _completed_job(settings, content: str):
    created = create_job(settings, JobCreateRequest(parent_smiles="CCO", comment="results"))
    repository = JobRepository(settings)
    assert repository.claim_next("test-worker") is not None
    path = settings.jobs_dir / str(created.id) / "output" / "202610161200_results.csv"
    path.write_text(content, encoding="utf-8")
    return repository.finish(
        created.id,
        target=JobStatus.COMPLETED,
        exit_code=0,
        result_csv_path=relative_to_repository(settings, path),
    )


def test_result_schema_filter_statistics_and_downloads(initialized_settings):
    settings = initialized_settings
    job = _completed_job(
        settings,
        "SMILES,CtRotBonds,NLL\nCCO,0,1.1\nCCCC,1,2.2\ninvalid,2,3.3\n",
    )
    with TestClient(create_app(settings)) as client:
        schema = client.get(f"/api/jobs/{job.id}/results/schema")
        molecules = client.get(
            f"/api/jobs/{job.id}/results/molecules?mw_min=50&sort_by=MW&sort_direction=desc"
        )
        statistics = client.get(f"/api/jobs/{job.id}/results/statistics?mw_min=50")
        csv_download = client.get(f"/api/jobs/{job.id}/results/download.csv?mw_min=50")
        sdf_download = client.get(f"/api/jobs/{job.id}/results/download.sdf?mw_min=50")

    assert schema.status_code == 200
    assert schema.json()["total_rows"] == 3
    assert schema.json()["valid_rows"] == 2
    assert schema.json()["invalid_smiles"] == 1
    assert {item["name"] for item in schema.json()["descriptors"]} >= {
        "MW", "RotB", "HBD", "HBA", "LogP", "TPSA", "QED", "SAScore",
        "FractionCSP3", "RingCount", "AromaticRingCount", "AliphaticRingCount",
        "HeavyAtomCount", "FormalCharge", "LipinskiViolations", "MolMR",
        "ParentSimilarity", "ScaffoldSimilarity", "NLL"
    }
    assert molecules.status_code == 200
    assert molecules.json()["total"] == 1
    assert molecules.json()["items"][0]["canonical_smiles"] == "CCCC"
    assert "<svg" in molecules.json()["items"][0]["structure_svg"]
    assert molecules.json()["items"][0]["values"]["QED"] is not None
    assert molecules.json()["items"][0]["values"]["ParentSimilarity"] is not None
    assert statistics.json()["filtered_rows"] == 1
    assert statistics.json()["summaries"]["MW"]["count"] == 1
    assert statistics.json()["scatter_x"] == "MW"
    assert statistics.json()["scatter_y"] == "LogP"
    assert len(statistics.json()["scatter"]) == 1
    assert csv_download.text.count("\n") == 2
    assert "QED" in csv_download.text.splitlines()[0]
    assert "ParentSimilarity" in csv_download.text.splitlines()[0]
    assert sdf_download.text.count("$$$$") == 1
    assert "202610161200_results.csv" in csv_download.headers["content-disposition"]


def test_results_require_completion_and_validate_query(initialized_settings):
    settings = initialized_settings
    completed = _completed_job(settings, "SMILES\nCCO\n")
    queued = create_job(settings, JobCreateRequest(parent_smiles="CCO"))
    with TestClient(create_app(settings)) as client:
        not_ready = client.get(f"/api/jobs/{queued.id}/results/schema")
        invalid_sort = client.get(
            f"/api/jobs/{completed.id}/results/molecules?sort_by=__class__"
        )
        invalid_range = client.get(
            f"/api/jobs/{completed.id}/results/molecules?mw_min=100&mw_max=10"
        )
        nonfinite_range = client.get(
            f"/api/jobs/{completed.id}/results/molecules?mw_min=nan"
        )
        advanced_filter = client.get(
            f"/api/jobs/{completed.id}/results/molecules?qed_min=0&qed_max=1"
        )

    assert not_ready.status_code == 409
    assert not_ready.json()["error"]["code"] == "results_not_ready"
    assert invalid_sort.status_code == 422
    assert invalid_range.status_code == 422
    assert nonfinite_range.status_code == 422
    assert advanced_filter.status_code == 200


def test_large_result_is_paginated_instead_of_returned_at_once(initialized_settings):
    rows = "".join(f"C,{index}\n" for index in range(10_001))
    job = _completed_job(initialized_settings, "SMILES,NLL\n" + rows)

    with TestClient(create_app(initialized_settings)) as client:
        response = client.get(f"/api/jobs/{job.id}/results/molecules?limit=50")

    assert response.status_code == 200
    assert response.json()["total"] == 10_001
    assert len(response.json()["items"]) == 50
