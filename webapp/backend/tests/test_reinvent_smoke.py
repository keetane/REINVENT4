from __future__ import annotations

import csv
import os
from dataclasses import replace
from pathlib import Path

import pytest

from app.db.jobs import JobRepository
from app.models.job import JobStatus
from app.schemas.jobs import JobCreateRequest
from app.services.jobs import create_job
from app.services.naming import confined_path
from worker.engine import Worker


@pytest.mark.reinvent
@pytest.mark.skipif(
    os.getenv("RUN_REINVENT_SMOKE") != "1",
    reason="set RUN_REINVENT_SMOKE=1 for the model-backed smoke test",
)
def test_real_libinvent_smoke(initialized_settings):
    source_root = Path(
        os.getenv("REINVENT_SMOKE_SOURCE", "/Users/keetane/Documents/modeling/REINVENT4")
    )
    runtime = Path(
        os.getenv(
            "REINVENT_SMOKE_PYTHON", "/Users/keetane/miniforge3/envs/r4/bin/python"
        )
    )
    settings = replace(
        initialized_settings,
        worker_runner="reinvent",
        reinvent_runtime_python=runtime,
        reinvent_source_root=source_root,
    )
    created = create_job(
        settings,
        JobCreateRequest(
            parent_smiles="N#CC[C@H](C1CCCC1)[*:1]",
            comment="libinvent smoke",
            run_conditions={
                "model_key": "libinvent",
                "num_smiles": 2,
                "device": "cpu",
                "unique_molecules": True,
                "randomize_smiles": True,
            },
        ),
    )

    Worker(settings, worker_id="real-smoke").run_once()

    completed = JobRepository(settings).get(created.id)
    assert completed.status == JobStatus.COMPLETED
    assert completed.exit_code == 0
    assert completed.model_path == "priors/libinvent.prior"
    assert completed.config_snapshot_path is not None
    assert completed.command_metadata["argv"][1:3] == ["-m", "reinvent"]
    result = confined_path(settings.repository_root, completed.result_csv_path)
    log = confined_path(settings.repository_root, completed.log_path)
    config = confined_path(settings.repository_root, completed.config_snapshot_path)
    with result.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert "SMILES" in rows[0]
    assert result.name.endswith("_libinvent_smoke.csv")
    assert "Finished REINVENT" in log.read_text(encoding="utf-8")
    assert config.is_file()
