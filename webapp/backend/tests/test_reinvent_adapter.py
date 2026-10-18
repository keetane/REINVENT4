from dataclasses import replace
from pathlib import Path

import pytest

from app.db.jobs import JobRepository
from app.errors import AppError
from app.reinvent_adapter import AdapterConfigurationError, build_reinvent_command
from app.schemas.jobs import JobCreateRequest
from app.services.jobs import create_job
from app.services.naming import confined_path
from worker.engine import Worker


def _adapter_settings(initialized_settings, tmp_path: Path):
    source_root = tmp_path / "source"
    (source_root / "reinvent").mkdir(parents=True)
    (source_root / "reinvent" / "__init__.py").write_text("", encoding="utf-8")
    (source_root / "priors").mkdir()
    (source_root / "priors" / "libinvent.prior").write_bytes(b"model")
    runtime = tmp_path / "runtime-python"
    runtime.write_text("#!/bin/sh\n", encoding="utf-8")
    runtime.chmod(0o755)
    return replace(
        initialized_settings,
        worker_runner="reinvent",
        reinvent_runtime_python=runtime,
        reinvent_source_root=source_root,
    )


def test_adapter_generates_allowlisted_argv_and_toml(initialized_settings, tmp_path):
    settings = _adapter_settings(initialized_settings, tmp_path)
    created = create_job(
        settings,
        JobCreateRequest(
            parent_smiles="N#CCC[*:1]",
            comment="hello; touch /tmp/nope",
            run_conditions={
                "model_key": "libinvent",
                "num_smiles": 2,
                "device": "cpu",
            },
        ),
    )
    job_dir = confined_path(settings.repository_root, created.job_dir)

    command = build_reinvent_command(settings, created, job_dir)

    assert command.argv[:3] == [str(settings.reinvent_runtime_python), "-m", "reinvent"]
    assert "hello" not in " ".join(command.argv)
    assert command.model_path == "priors/libinvent.prior"
    toml = command.config_path.read_text(encoding="utf-8")
    assert 'run_type = "sampling"' in toml
    assert "num_smiles = 2" in toml
    assert "touch /tmp/nope" not in toml


def test_adapter_rejects_client_selected_model_path(initialized_settings, tmp_path):
    settings = _adapter_settings(initialized_settings, tmp_path)
    with pytest.raises(AppError, match="利用できません"):
        create_job(
            settings,
            JobCreateRequest(
                parent_smiles="CCO",
                run_conditions={"model_key": "../../private/model"},
            ),
        )


def test_worker_records_adapter_configuration_failure(initialized_settings, tmp_path):
    settings = replace(
        initialized_settings,
        worker_runner="reinvent",
        reinvent_runtime_python=tmp_path / "missing-python",
        reinvent_source_root=tmp_path / "missing-source",
    )
    created = create_job(
        settings,
        JobCreateRequest(parent_smiles="CCO[*]", run_conditions={"model_key": "libinvent"}),
    )

    Worker(settings).run_once()

    failed = JobRepository(settings).get(created.id)
    assert failed.status.value == "failed"
    assert failed.error_kind == "adapter_configuration"
    assert "Python is not executable" in failed.error_summary


def test_adapter_rejects_non_boolean_sampling_flags(initialized_settings, tmp_path):
    settings = _adapter_settings(initialized_settings, tmp_path)
    created = create_job(
        settings,
        JobCreateRequest(
            parent_smiles="CCO[*]",
            run_conditions={
                "model_key": "libinvent",
                "unique_molecules": "false",
            },
        ),
    )
    with pytest.raises(AdapterConfigurationError, match="must be a boolean"):
        build_reinvent_command(
            settings,
            created,
            confined_path(settings.repository_root, created.job_dir),
        )
