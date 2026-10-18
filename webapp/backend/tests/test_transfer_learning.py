from __future__ import annotations

from dataclasses import replace

import pytest
from fastapi.testclient import TestClient

from app.db.jobs import JobRepository
from app.db.models import ModelRepository
from app.constants import MAX_DATA_UPLOAD_BYTES
from app.main import create_app
from app.reinvent_adapter import build_transfer_learning_command
from app.schemas.training import TrainingJobCreateRequest
from app.services.naming import confined_path
from app.services.training import create_training_job
from worker.engine import Worker


IBUPROFEN = "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O"
NAPROXEN = "COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1"
LIDOCAINE = "CCN(CC)C(=O)c1c(C)cccc1C"
WARFARIN = "CC(=O)CC(c1ccccc1)c1c(O)oc2ccccc12"
CSV_DATA = f"SMILES,name\n{IBUPROFEN},ibuprofen\n{NAPROXEN},naproxen\n{IBUPROFEN},duplicate\nnot-smiles,bad\n{LIDOCAINE},lidocaine\n"
CHEMBL_CSV_DATA = f'Smiles;"pChEMBL Value";Molecule Name\n{IBUPROFEN};5.1;ibuprofen\n{NAPROXEN};6.2;naproxen\n'


def test_training_data_inspection_filters_invalid_and_duplicate_smiles(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/training-data/inspect",
            json={"file_name": "training.csv", "content": CSV_DATA, "smiles_column": "SMILES"},
        )

    assert response.status_code == 200
    assert response.json()["columns"] == ["SMILES", "name"]
    assert response.json()["total_rows"] == 5
    assert response.json()["valid_count"] == 3
    assert response.json()["invalid_count"] == 1
    assert response.json()["duplicate_count"] == 1
    assert response.json()["standardization_version"] == "1.0.0"


def test_training_data_accepts_chembl_semicolon_csv_and_larger_payloads(settings):
    assert MAX_DATA_UPLOAD_BYTES == 52_428_800
    request = TrainingJobCreateRequest(
        file_name="chembl.csv",
        content=CHEMBL_CSV_DATA + (" " * 10_000_000),
        smiles_column="Smiles",
        base_model="reinvent",
        comment="chembl",
    )

    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/training-data/inspect",
            json={"file_name": "chembl.csv", "content": CHEMBL_CSV_DATA, "smiles_column": "Smiles"},
        )

    assert request.smiles_column == "Smiles"
    assert response.status_code == 200
    assert response.json()["columns"] == ["Smiles", "pChEMBL Value", "Molecule Name"]
    assert response.json()["valid_count"] == 2


def test_training_job_completes_and_registers_persistent_model(settings):
    with TestClient(create_app(settings)) as client:
        created = client.post(
            "/api/training-jobs",
            json={
                "file_name": "training.csv",
                "content": CSV_DATA,
                "smiles_column": "SMILES",
                "base_model": "reinvent",
                "comment": "kinase_A",
                "train_fraction": 0.8,
                "num_epochs": 2,
                "batch_size": 2,
                "learning_rate": 0.0001,
                "device": "cpu",
            },
        )
        assert created.status_code == 202
        assert created.json()["job_type"] == "transfer_learning"
        job_id = created.json()["id"]

    assert Worker(settings, worker_id="tl-worker").run_once()
    completed = JobRepository(settings).get(job_id)
    assert completed.status.value == "completed"
    models = ModelRepository(settings).list()
    assert len(models) == 1
    assert models[0].name.endswith("_kinase_A.model")
    assert (settings.repository_root / models[0].file_path).read_bytes() == b"fake trained model"
    assert (settings.repository_root / models[0].config_path).is_file()
    assert (settings.repository_root / models[0].log_path).is_file()

    with TestClient(create_app(settings)) as restarted_client:
        listed = restarted_client.get("/api/models")
    custom = [item for item in listed.json()["items"] if item["source"] == "transfer_learning"]
    assert len(custom) == 1
    assert custom[0]["training_data_count"] == 2
    assert custom[0]["training_conditions"]["smiles_standardization"]["version"] == "1.0.0"
    with TestClient(create_app(settings)) as client:
        assert client.get(f"/api/jobs/{job_id}/standardization/accepted").status_code == 200
        assert client.get(f"/api/jobs/{job_id}/standardization/rejected").status_code == 200


def test_training_progress_preserves_one_history_point_per_epoch(initialized_settings):
    settings = initialized_settings
    job = create_training_job(
        settings,
        TrainingJobCreateRequest(
            file_name="training.smi",
            content=f"{IBUPROFEN}\n{NAPROXEN}\n{LIDOCAINE}\n",
            base_model="reinvent",
            comment="history",
            num_epochs=3,
            batch_size=2,
            learning_rate=0.0001,
            device="cpu",
        ),
    )
    repository = JobRepository(settings)
    assert repository.claim_next("test-worker") is not None

    repository.set_training_progress(
        job.id, epoch=1, total_epochs=3, training_loss=1.0, validation_loss=1.1
    )
    repository.set_training_progress(
        job.id, epoch=1, total_epochs=3, training_loss=0.9, validation_loss=1.05
    )
    updated = repository.set_training_progress(
        job.id, epoch=2, total_epochs=3, training_loss=0.7, validation_loss=0.8
    )

    assert updated.metrics["history"] == [
        {"epoch": 1, "training_loss": 0.9, "validation_loss": 1.05},
        {"epoch": 2, "training_loss": 0.7, "validation_loss": 0.8},
    ]


def test_reinvent_generation_does_not_require_parent_and_snapshots_model(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/jobs",
            json={
                "parent_smiles": "",
                "comment": "prior",
                "run_conditions": {"model_id": "standard:reinvent"},
            },
        )

    assert response.status_code == 202
    job = JobRepository(settings).get(response.json()["id"])
    assert job.parent_smiles_input == ""
    assert job.run_conditions["resolved_model_path"] == "priors/reinvent.prior"


def test_failed_training_does_not_register_model(initialized_settings, tmp_path):
    settings = replace(
        initialized_settings,
        worker_runner="reinvent",
        reinvent_runtime_python=tmp_path / "missing-python",
        reinvent_source_root=tmp_path / "missing-source",
    )
    job = create_training_job(
        settings,
        TrainingJobCreateRequest(
            file_name="training.smi",
            content=f"{IBUPROFEN}\n{NAPROXEN}\n{LIDOCAINE}\n",
            base_model="reinvent",
            comment="failure",
            num_epochs=1,
            batch_size=2,
            learning_rate=0.0001,
            device="cpu",
        ),
    )

    Worker(settings).run_once()

    failed = JobRepository(settings).get(job.id)
    assert failed.status.value == "failed"
    assert ModelRepository(settings).list() == []
    assert not (settings.repository_root / failed.run_conditions["output_model_path"]).exists()


def test_model_in_use_cannot_be_deleted(settings):
    with TestClient(create_app(settings)) as client:
        created = client.post(
            "/api/training-jobs",
            json={
                "file_name": "training.smi",
                "content": f"{IBUPROFEN}\n{NAPROXEN}\n{LIDOCAINE}\n",
                "base_model": "reinvent",
                "comment": "delete_guard",
                "num_epochs": 1,
                "batch_size": 2,
                "learning_rate": 0.0001,
                "device": "cpu",
            },
        )
        Worker(settings).run_once()
        model = ModelRepository(settings).list()[0]
        generation = client.post(
            "/api/jobs",
            json={
                "parent_smiles": "",
                "run_conditions": {"model_id": str(model.id)},
            },
        )
        assert generation.status_code == 202
        deleted = client.delete(f"/api/models/{model.id}")

    assert deleted.status_code == 409
    assert deleted.json()["error"]["code"] == "model_in_use"


@pytest.mark.parametrize(
    ("base_model", "prior_name", "expected_pairs"),
    [
        ("reinvent", "reinvent.prior", False),
        ("mol2mol", "mol2mol_scaffold_generic.prior", True),
    ],
)
def test_transfer_learning_toml_matches_supported_reinvent_parameters(
    initialized_settings, tmp_path, base_model, prior_name, expected_pairs
):
    source_root = tmp_path / "source"
    (source_root / "reinvent").mkdir(parents=True)
    (source_root / "reinvent" / "__init__.py").write_text("", encoding="utf-8")
    (source_root / "priors").mkdir()
    (source_root / "priors" / prior_name).write_bytes(b"prior")
    runtime = tmp_path / "python"
    runtime.write_text("#!/bin/sh\n", encoding="utf-8")
    runtime.chmod(0o755)
    settings = replace(
        initialized_settings,
        reinvent_runtime_python=runtime,
        reinvent_source_root=source_root,
    )
    job = create_training_job(
        settings,
        TrainingJobCreateRequest(
            file_name="training.smi",
            content=f"{IBUPROFEN}\n{NAPROXEN}\n{LIDOCAINE}\n{WARFARIN}\n",
            base_model=base_model,
            comment=f"{base_model}_test",
            num_epochs=2,
            batch_size=2,
            learning_rate=0.0002,
            tanimoto_lower=0.6,
            tanimoto_upper=0.9,
            device="cpu",
        ),
    )
    command = build_transfer_learning_command(
        settings, job, confined_path(settings.repository_root, job.job_dir)
    )
    toml = command.config_path.read_text(encoding="utf-8")
    assert 'run_type = "transfer_learning"' in toml
    assert "num_epochs = 2" in toml
    assert "batch_size = 2" in toml
    assert "lr = 0.0002" in toml
    assert "[responder]" in toml
    assert ('pairs.type = "tanimoto"' in toml) is expected_pairs
