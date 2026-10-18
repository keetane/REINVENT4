from __future__ import annotations

from dataclasses import replace

from fastapi.testclient import TestClient

from app.db.jobs import JobRepository
from app.db.models import ModelRepository
from app.main import create_app
from app.reinvent_adapter import build_rl_command
from app.schemas.rl import RlJobCreateRequest, RlPropertyComponent
from app.services.naming import confined_path
from app.services.rl import create_rl_job
from worker.engine import Worker


def rl_payload(**updates):
    value = {
        "generation_model_id": "standard:reinvent",
        "comment": "qed_rl",
        "device": "cpu",
        "property_components": [{"component": "QED", "weight": 1, "transform": "none"}],
        "batch_size": 8,
        "sigma": 128,
        "learning_rate": 0.0001,
        "minimum_steps": 1,
        "maximum_steps": 2,
        "target_mean_score": 0.8,
        "diversity_filter": True,
    }
    value.update(updates)
    return value


def test_rl_job_completes_registers_model_and_can_generate(settings):
    with TestClient(create_app(settings)) as client:
        created = client.post("/api/rl-jobs", json=rl_payload())
        assert created.status_code == 202
        assert created.json()["job_type"] == "reinforcement_learning"
        job_id = created.json()["id"]

    assert Worker(settings, worker_id="rl-worker").run_once()
    completed = JobRepository(settings).get(job_id)
    assert completed.status.value == "completed"
    assert completed.result_csv_path is not None
    audit = (settings.repository_root / completed.result_csv_path).read_text(encoding="utf-8")
    assert "Generated_SMILES,Standardized_SMILES,MolecularWeight,FormalCharge" in audit
    models = ModelRepository(settings).list()
    assert len(models) == 1
    assert models[0].training_conditions["source"] == "reinforcement_learning"

    with TestClient(create_app(settings)) as client:
        listed = client.get("/api/models")
        generation = client.post(
            "/api/jobs",
            json={"parent_smiles": "", "comment": "from_rl", "run_conditions": {"model_id": str(models[0].id)}},
        )
    custom = [item for item in listed.json()["items"] if item["id"] == str(models[0].id)][0]
    assert custom["source"] == "reinforcement_learning"
    assert generation.status_code == 202


def test_rl_toml_uses_official_staged_learning_fields(initialized_settings, tmp_path):
    source = tmp_path / "source"
    (source / "reinvent").mkdir(parents=True)
    (source / "reinvent" / "__init__.py").write_text("", encoding="utf-8")
    (source / "priors").mkdir()
    (source / "priors" / "reinvent.prior").write_bytes(b"prior")
    runtime = tmp_path / "python"
    runtime.write_text("#!/bin/sh\n", encoding="utf-8")
    runtime.chmod(0o755)
    settings = replace(initialized_settings, reinvent_source_root=source, reinvent_runtime_python=runtime)
    job = create_rl_job(
        settings,
        RlJobCreateRequest(
            generation_model_id="standard:reinvent",
            comment="config",
            device="cpu",
            property_components=[RlPropertyComponent(component="QED")],
            minimum_steps=1,
            maximum_steps=2,
        ),
    )
    command = build_rl_command(settings, job, confined_path(settings.repository_root, job.job_dir))
    toml = command.config_path.read_text(encoding="utf-8")
    assert 'run_type = "staged_learning"' in toml
    assert '[learning_strategy]' in toml
    assert 'type = "dap"' in toml
    assert '[[stage]]' in toml
    assert '[stage.scoring.component.QED]' in toml
    assert '[responder]' in toml


def test_mol2mol_rl_requires_input_smiles(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/rl-jobs",
            json=rl_payload(generation_model_id="standard:mol2mol_scaffold_generic"),
        )
    assert response.status_code == 422
    assert response.json()["error"]["code"] == "rl_input_required"


def test_rl_preview_includes_properties_and_total_score(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/rl/preview",
            json={
                "smiles": ["CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "CCO", "not-smiles"],
                "property_components": [
                    {"component": "QED", "weight": 2, "transform": "none"},
                    {"component": "MolecularWeight", "weight": 1, "transform": "double_sigmoid", "low": 200, "high": 500},
                ],
                "aggregation": "geometric_mean",
            },
        )
    assert response.status_code == 200
    rows = response.json()["rows"]
    assert rows[0]["valid"] is True
    assert [item["kind"] for item in rows[0]["components"]] == ["property", "property"]
    assert 0 < rows[0]["total_score"] <= 1
    assert rows[1]["valid"] is False
    assert rows[1]["exclusion_reason"] == "molecular_weight_not_above_200"
    assert rows[1]["total_score"] == 0
    assert rows[2]["valid"] is False
    assert rows[2]["total_score"] == 0
