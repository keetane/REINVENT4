from types import SimpleNamespace
from uuid import uuid4

import pytest

from worker.__main__ import initialize_worker
from worker.engine import Worker


def test_worker_can_initialize_independently(settings):
    initialize_worker(settings)

    assert settings.database_path.is_file()
    assert settings.jobs_dir.is_dir()


def test_chemprop_checkpoint_name_is_not_parsed_as_validation_loss(settings, tmp_path):
    updates = []
    worker = Worker(settings)
    worker.repository = SimpleNamespace(
        set_training_progress=lambda *args, **kwargs: updates.append((args, kwargs))
    )
    log_path = tmp_path / "chemprop.log"
    log_path.write_text(
        "Epoch 24/29 train_loss_epoch: 0.421 val_loss: 0.468\n"
        "Restoring best-epoch=14-val_loss=0.42.ckpt\n",
        encoding="utf-8",
    )
    job = SimpleNamespace(
        id=uuid4(),
        job_type="qsar_training",
        run_conditions={"algorithm": "chemprop_v2", "hyperparameters": {"epochs": 30}},
    )

    assert worker._update_log_progress(job, log_path, -1) == 25
    assert updates[0][1] == {
        "epoch": 25,
        "total_epochs": 30,
        "training_loss": pytest.approx(0.421),
        "validation_loss": pytest.approx(0.468),
    }
