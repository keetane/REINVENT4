from __future__ import annotations

import csv
from io import StringIO

import pytest
from fastapi.testclient import TestClient

from app.db.jobs import JobRepository
from app.db.qsar_models import QsarModelRepository
from app.errors import AppError
from app.main import create_app
from app.schemas.qsar import QsarDataRequest
from app.services.qsar import process_qsar_data, validate_qsar_device
from worker.engine import Worker


CSV_DATA = """SMILES,Activity
CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O,100
COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1,200
CCN(CC)C(=O)c1c(C)cccc1C,300
CC(=O)CC(c1ccccc1)c1c(O)oc2ccccc12,400
CN(C)CCOC(c1ccccc1)c1ccccc1,500
CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1,600
CC1=NN(C2=CC=C(C=C2)S(N)(=O)=O)C(C)=C1C3=CC=CC=C3,700
CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O,900
not-smiles,50
CCOC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+]([O-])=O,
"""
CHEMBL_CSV_DATA = '''Smiles;"pChEMBL Value";Molecule Name
CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O;5.1;ibuprofen
COc1ccc2cc([C@@H](C)C(=O)O)ccc2c1;6.2;naproxen
CCN(CC)C(=O)c1c(C)cccc1C;;lidocaine
'''


def payload(**updates):
    value = {
        "file_name": "activity.csv",
        "content": CSV_DATA,
        "smiles_column": "SMILES",
        "target_column": "Activity",
        "task_type": "regression",
        "target_transform": "p_activity",
        "concentration_unit": "nM",
        "duplicate_strategy": "exclude_conflicts",
    }
    value.update(updates)
    return value


def test_qsar_preprocessing_validates_converts_and_reports_conflicts(settings):
    processed = process_qsar_data(QsarDataRequest(**payload()), max_atoms=200)

    assert processed.total_rows == 10
    assert len(processed.rows) == 6
    assert processed.invalid_smiles_count == 1
    assert processed.missing_target_count == 1
    assert processed.duplicate_count == 1
    assert processed.conflicting_duplicate_count == 2
    assert {row.reason for row in processed.rejected_rows} >= {
        "invalid_smiles", "missing_target", "conflicting_duplicate"
    }
    assert next(target for smiles, target in processed.rows if "COc1ccc2" in smiles) == pytest.approx(6.69897)


def test_qsar_preprocessing_accepts_chembl_semicolon_csv():
    processed = process_qsar_data(
        QsarDataRequest(
            file_name="chembl.csv",
            content=CHEMBL_CSV_DATA,
            smiles_column="Smiles",
            target_column="pChEMBL Value",
        ),
        max_atoms=200,
    )

    assert processed.columns == ["Smiles", "pChEMBL Value", "Molecule Name"]
    assert processed.total_rows == 3
    assert len(processed.rows) == 2
    assert [target for _, target in processed.rows] == [5.1, 6.2]
    assert processed.missing_target_count == 1


def test_qsar_training_job_registers_persistent_artifacts(settings):
    with TestClient(create_app(settings)) as client:
        inspected = client.post("/api/qsar-data/inspect", json=payload(duplicate_strategy="keep_first"))
        assert inspected.status_code == 200
        assert inspected.json()["valid_count"] == 7
        created = client.post(
            "/api/qsar-jobs",
            json={
                **payload(duplicate_strategy="keep_first"),
                "algorithm": "random_forest",
                "comment": "activity_rf",
                "fingerprint": "morgan",
                "feature_mode": "fingerprint",
                "split_method": "random",
                "train_fraction": 0.8,
                "cross_validation_folds": 3,
                "ad_threshold": 0.4,
                "device": "cpu",
            },
        )
        assert created.status_code == 202
        job_id = created.json()["id"]
        assert created.json()["job_type"] == "qsar_training"

    assert Worker(settings, worker_id="qsar-worker").run_once()
    job = JobRepository(settings).get(job_id)
    assert job.status.value == "completed"
    assert job.metrics["test"]["r2"] == 0.7
    models = QsarModelRepository(settings).list()
    assert len(models) == 1
    model = models[0]
    assert model.algorithm == "random_forest"
    assert model.metadata["smiles_standardization"]["version"] == "1.0.0"
    assert (settings.repository_root / model.file_path).read_bytes() == b"trusted fake QSAR model"

    with TestClient(create_app(settings)) as client:
        listed = client.get("/api/qsar-models")
        evaluation = client.get(f"/api/qsar-models/{model.id}/artifacts/evaluation")
        deleted = client.delete(f"/api/qsar-models/{model.id}")
        accepted = client.get(f"/api/jobs/{job_id}/standardization/accepted")
        rejected = client.get(f"/api/jobs/{job_id}/standardization/rejected")
    assert listed.status_code == 200
    assert listed.json()["items"][0]["name"].endswith("_activity_rf")
    assert evaluation.status_code == 200
    assert accepted.status_code == 200
    assert rejected.status_code == 200
    assert deleted.status_code == 200
    assert not (settings.repository_root / model.file_path).exists()


def test_qsar_model_upload_infers_algorithm_and_comment(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/qsar-models/upload",
            content=b"uploaded-joblib",
            headers={
                "Content-Type": "application/octet-stream",
                "X-Model-Filename": "20261019_QSAR_RandomForest_chembl_p2x3.joblib",
            },
        )

    assert response.status_code == 201
    model = response.json()
    assert model["name"] == "chembl_p2x3"
    assert model["algorithm"] == "random_forest"
    assert model["metadata"]["uploaded"] is True


def test_unavailable_qsar_algorithms_are_rejected(settings):
    with TestClient(create_app(settings)) as client:
        response = client.post(
            "/api/qsar-jobs",
            json={**payload(duplicate_strategy="keep_first"), "algorithm": "chemprop_v2", "comment": "chemprop"},
        )
    assert response.status_code == 422


def test_qsar_gpu_device_is_limited_to_native_gpu_algorithms():
    validate_qsar_device("cuda", ["lightgbm", "xgboost", "catboost"])
    with pytest.raises(AppError, match="CUDA"):
        validate_qsar_device("cuda", ["mlp"])
    with pytest.raises(AppError, match="MPS"):
        validate_qsar_device("mps", ["random_forest"])
