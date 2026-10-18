from __future__ import annotations

import os
from pathlib import Path

import pytest

from app.config import Settings
from app.db.migrations import apply_migrations
from app.db.jobs import JobRepository
from app.db.qsar_models import QsarModelRepository
from app.schemas.qsar import QsarJobCreateRequest
from app.services.qsar import create_qsar_job
from worker.engine import Worker


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
CHEMPROP_PYTHON = REPOSITORY_ROOT / ".runtimes" / "chemprop-v2" / "bin" / "python"
QSARTUNA_PYTHON = REPOSITORY_ROOT / ".runtimes" / "qsartuna" / "bin" / "python"
REINVENT_PYTHON = Path(os.getenv("REINVENT_RUNTIME_PYTHON", REPOSITORY_ROOT / ".runtimes" / "reinvent" / "bin" / "python"))
if not REINVENT_PYTHON.is_file():
    REINVENT_PYTHON = Path.home() / "miniforge3" / "envs" / "r4" / "bin" / "python"


def _training_csv() -> str:
    rows = ["SMILES,Activity"]
    rows.extend(f"{'C' * carbon_count},{5 + carbon_count / 20:.3f}" for carbon_count in range(15, 45))
    return "\n".join(rows) + "\n"


@pytest.mark.external_qsar
@pytest.mark.parametrize(
    ("algorithm", "hyperparameters"),
    [
        ("chemprop_v2", {"epochs": 1, "batch_size": 8, "learning_rate": 0.0001, "patience": 1}),
        ("qsartuna", {"n_trials": 1, "n_startup_trials": 0, "timeout_seconds": 180}),
        ("mlp", {"hidden_layer_sizes": [16], "alpha": 0.0001, "learning_rate_init": 0.001, "max_iter": 20}),
    ],
)
def test_real_qsar_algorithm_completes_and_registers(tmp_path: Path, algorithm: str, hyperparameters: dict):
    if os.getenv("REINVENT_QSAR_EXTERNAL_SMOKE") != "1":
        pytest.skip("set REINVENT_QSAR_EXTERNAL_SMOKE=1 to run real QSAR runtimes")
    if algorithm == "chemprop_v2" and not CHEMPROP_PYTHON.is_file():
        pytest.skip("Chemprop runtime is unavailable")
    if algorithm == "qsartuna" and not QSARTUNA_PYTHON.is_file():
        pytest.skip("QSARtuna runtime is unavailable")
    if not REINVENT_PYTHON.is_file():
        pytest.skip("REINVENT runtime is unavailable")

    settings = Settings(
        repository_root=tmp_path,
        database_path=tmp_path / "jobs.db",
        jobs_dir=tmp_path / "jobs",
        models_dir=tmp_path / "models",
        worker_runner="reinvent",
        reinvent_runtime_python=REINVENT_PYTHON,
        reinvent_source_root=REPOSITORY_ROOT,
        chemprop_runtime_python=CHEMPROP_PYTHON,
        qsartuna_runtime_python=QSARTUNA_PYTHON,
        worker_poll_interval_seconds=0.01,
    )
    apply_migrations(settings)
    settings.jobs_dir.mkdir(parents=True, exist_ok=True)
    job = create_qsar_job(settings, QsarJobCreateRequest(
        file_name="external-smoke.csv",
        content=_training_csv(),
        smiles_column="SMILES",
        target_column="Activity",
        task_type="regression",
        algorithm=algorithm,
        comment=f"smoke_{algorithm}",
        split_method="random",
        train_fraction=0.8,
        cross_validation_folds=2,
        bit_length=128,
        device="cpu",
        advanced=True,
        hyperparameters=hyperparameters,
    ))

    assert Worker(settings, worker_id=f"smoke-{algorithm}").run_once()
    completed = JobRepository(settings).get(job.id)
    log_text = (settings.repository_root / completed.log_path).read_text(encoding="utf-8", errors="replace") if completed.log_path else ""
    assert completed.status.value == "completed", f"{completed.error_kind}: {completed.error_summary}\n{log_text}"
    assert completed.exit_code == 0
    models = QsarModelRepository(settings).list()
    assert len(models) == 1
    assert models[0].algorithm == algorithm
    assert (settings.repository_root / models[0].file_path).stat().st_size > 0
    assert (settings.repository_root / models[0].evaluation_path).is_file()
    assert (settings.repository_root / models[0].predictions_path).is_file()
