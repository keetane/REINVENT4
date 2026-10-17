from __future__ import annotations

import csv
import threading
import time
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app.db.jobs import JobRepository
from app.main import create_app
from app.models.job import JobStatus
from app.schemas.jobs import JobCreateRequest
from app.services.jobs import create_job
from app.services.naming import confined_path
from worker.engine import Worker


def _wait_for_status(repository, job_id, expected, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        job = repository.get(job_id)
        if job.status == expected:
            return job
        time.sleep(0.02)
    raise AssertionError(f"job did not reach {expected}")


def test_worker_completes_after_api_client_disconnect(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/jobs",
            json={
                "parent_smiles": "CCO",
                "comment": "disconnect",
                "run_conditions": {"fake_rows": 3},
            },
        )
        job_id = response.json()["id"]

    assert Worker(settings, worker_id="integration-worker").run_once()

    job = JobRepository(settings).get(job_id)
    assert job.status == JobStatus.COMPLETED
    assert job.exit_code == 0
    result = confined_path(settings.repository_root, job.result_csv_path)
    with result.open(newline="", encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle))) == 3
    assert result.name.endswith("_disconnect.csv")


def test_api_restart_does_not_interrupt_independent_worker(settings):
    with TestClient(create_app(settings)) as first_client:
        created = first_client.post(
            "/api/jobs",
            json={
                "parent_smiles": "CCO",
                "comment": "restart",
                "run_conditions": {"fake_delay_seconds": 0.4},
            },
        ).json()

    thread = threading.Thread(target=Worker(settings).run_once)
    thread.start()
    repository = JobRepository(settings)
    _wait_for_status(repository, created["id"], JobStatus.RUNNING)

    with TestClient(create_app(settings)) as restarted_client:
        running = restarted_client.get(f"/api/jobs/{created['id']}")
        assert running.status_code == 200
        assert running.json()["status"] == "running"

    thread.join(timeout=5)
    assert not thread.is_alive()
    assert repository.get(created["id"]).status == JobStatus.COMPLETED


def test_cancel_stops_only_the_running_job(initialized_settings):
    slow = create_job(
        initialized_settings,
        JobCreateRequest(
            parent_smiles="CCO",
            comment="slow",
            run_conditions={"fake_delay_seconds": 5},
        ),
    )
    queued = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCC", comment="queued"),
    )
    repository = JobRepository(initialized_settings)
    thread = threading.Thread(target=Worker(initialized_settings).run_once)
    thread.start()
    _wait_for_status(repository, slow.id, JobStatus.RUNNING)

    repository.request_cancel(slow.id)
    thread.join(timeout=5)

    assert not thread.is_alive()
    assert repository.get(slow.id).status == JobStatus.CANCELLED
    assert repository.get(queued.id).status == JobStatus.QUEUED


def test_nonzero_exit_is_failed_with_original_diagnostic(initialized_settings):
    created = create_job(
        initialized_settings,
        JobCreateRequest(
            parent_smiles="CCO",
            comment="failure",
            run_conditions={"fake_exit_code": 7},
        ),
    )

    Worker(initialized_settings).run_once()

    job = JobRepository(initialized_settings).get(created.id)
    assert job.status == JobStatus.FAILED
    assert job.exit_code == 7
    assert job.error_kind == "nonzero_exit"
    assert "code 7" in job.error_summary


def test_same_minute_and_comment_never_overwrite_another_job(initialized_settings):
    first = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCO", comment="same"),
    )
    second = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCC", comment="same"),
    )
    worker = Worker(initialized_settings)
    assert worker.run_once()
    assert worker.run_once()

    repository = JobRepository(initialized_settings)
    first_done = repository.get(first.id)
    second_done = repository.get(second.id)
    assert first_done.result_csv_path != second_done.result_csv_path
    assert confined_path(
        initialized_settings.repository_root, first_done.result_csv_path
    ).is_file()
    assert confined_path(
        initialized_settings.repository_root, second_done.result_csv_path
    ).is_file()


@pytest.mark.parametrize(
    ("mode", "message"),
    [
        ("missing", "missing"),
        ("empty", "empty"),
        ("missing-smiles", "SMILES"),
    ],
)
def test_missing_or_invalid_worker_output_is_failed(
    initialized_settings, mode, message
):
    created = create_job(
        initialized_settings,
        JobCreateRequest(
            parent_smiles="CCO",
            comment=mode,
            run_conditions={"fake_output_mode": mode},
        ),
    )

    Worker(initialized_settings).run_once()

    failed = JobRepository(initialized_settings).get(created.id)
    assert failed.status == JobStatus.FAILED
    assert failed.error_kind == "invalid_output"
    assert message in failed.error_summary


def test_bounded_log_endpoint(initialized_settings):
    created = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCO", comment="logs"),
    )
    Worker(initialized_settings).run_once()

    with TestClient(create_app(initialized_settings)) as client:
        response = client.get(f"/api/jobs/{created.id}/logs?tail=1")

    assert response.status_code == 200
    assert response.json()["lines"] == ["fake worker finished"]
    assert response.json()["truncated"] is True
