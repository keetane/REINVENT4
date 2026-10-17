from pathlib import Path

from fastapi.testclient import TestClient

from app.main import create_app


def test_create_list_get_and_cancel_queued_job(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/jobs",
            json={
                "parent_smiles": "CCO",
                "comment": "../日本語 run",
                "run_conditions": {"fake_rows": 3},
            },
        )
        assert response.status_code == 202
        created = response.json()
        assert created["status"] == "queued"

        listed = client.get("/api/jobs?limit=10&offset=0").json()
        assert listed["total"] == 1
        assert listed["items"][0]["id"] == created["id"]

        detail = client.get(f"/api/jobs/{created['id']}")
        assert detail.status_code == 200

        cancelled = client.post(f"/api/jobs/{created['id']}/cancel")
        assert cancelled.status_code == 200
        assert cancelled.json()["status"] == "cancelled"

    workspace = Path(settings.jobs_dir, created["id"])
    assert workspace.is_dir()
    assert (workspace / "input" / "parent.smi").read_text() == "CCO\n"
    assert (workspace / "input" / "request.json").is_file()


def test_validation_and_not_found_errors_are_uniform(settings):
    with TestClient(create_app(settings)) as client:
        invalid = client.post("/api/jobs", json={"parent_smiles": ""})
        assert invalid.status_code == 422
        assert invalid.json()["error"]["code"] == "validation_error"

        missing = client.get("/api/jobs/12345678-1234-5678-1234-567812345678")
        assert missing.status_code == 404
        assert missing.json()["error"]["code"] == "job_not_found"

