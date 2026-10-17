from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace

from app.db.jobs import JobRepository
from app.schemas.jobs import JobCreateRequest
from app.services.jobs import create_job
from worker.engine import reconcile_active_jobs


def test_two_workers_cannot_claim_the_same_job(initialized_settings):
    created = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCO", comment="claim"),
    )
    repository = JobRepository(initialized_settings)

    with ThreadPoolExecutor(max_workers=2) as executor:
        claimed = list(
            executor.map(repository.claim_next, ["worker-one", "worker-two"])
        )

    actual = [job for job in claimed if job is not None]
    assert len(actual) == 1
    assert actual[0].id == created.id


def test_startup_reconciliation_interrupts_unowned_running_job(initialized_settings):
    created = create_job(
        initialized_settings,
        JobCreateRequest(parent_smiles="CCO", comment="reconcile"),
    )
    repository = JobRepository(initialized_settings)
    claimed = repository.claim_next("vanished-worker")
    assert claimed is not None

    interrupted = reconcile_active_jobs(initialized_settings)

    assert interrupted == [str(created.id)]
    assert repository.get(created.id).status.value == "interrupted"


def test_reconciliation_grace_avoids_claim_startup_race(initialized_settings):
    settings = replace(initialized_settings, reconciliation_grace_seconds=30.0)
    created = create_job(
        settings,
        JobCreateRequest(parent_smiles="CCO", comment="claim-race"),
    )
    repository = JobRepository(settings)
    assert repository.claim_next("starting-worker") is not None

    assert reconcile_active_jobs(settings) == []
    assert repository.get(created.id).status.value == "running"
