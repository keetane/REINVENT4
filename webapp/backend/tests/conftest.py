from pathlib import Path

import pytest

from app.config import Settings
from app.db.migrations import apply_migrations


@pytest.fixture
def settings(tmp_path: Path) -> Settings:
    return Settings(
        repository_root=tmp_path,
        database_path=tmp_path / "jobs.db",
        jobs_dir=tmp_path / "jobs",
        worker_poll_interval_seconds=0.01,
        cancel_grace_seconds=0.1,
        reconciliation_grace_seconds=0.0,
    )


@pytest.fixture
def initialized_settings(settings: Settings) -> Settings:
    apply_migrations(settings)
    settings.jobs_dir.mkdir(parents=True, exist_ok=True)
    return settings
