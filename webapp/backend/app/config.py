from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def _resolve_path(value: str | None, default: Path) -> Path:
    return Path(value).expanduser().resolve() if value else default.resolve()


@dataclass(frozen=True, slots=True)
class Settings:
    repository_root: Path
    database_path: Path
    jobs_dir: Path
    app_timezone: str = "Asia/Tokyo"
    cors_origins: tuple[str, ...] = ("http://127.0.0.1:5173", "http://localhost:5173")
    sqlite_busy_timeout_ms: int = 5_000
    worker_poll_interval_seconds: float = 2.0
    cancel_grace_seconds: float = 5.0
    reconciliation_grace_seconds: float = 30.0
    parent_max_atoms: int = 200
    worker_runner: str = "fake"
    reinvent_runtime_python: Path = Path(sys.executable)
    reinvent_source_root: Path = REPOSITORY_ROOT

    @classmethod
    def from_env(cls) -> "Settings":
        repository_root = _resolve_path(
            os.getenv("REINVENT_WEB_REPOSITORY_ROOT"), REPOSITORY_ROOT
        )
        origins = tuple(
            origin.strip()
            for origin in os.getenv(
                "REINVENT_WEB_CORS_ORIGINS",
                "http://127.0.0.1:5173,http://localhost:5173",
            ).split(",")
            if origin.strip()
        )
        return cls(
            repository_root=repository_root,
            database_path=_resolve_path(
                os.getenv("REINVENT_WEB_DATABASE_PATH"), repository_root / "jobs.db"
            ),
            jobs_dir=_resolve_path(
                os.getenv("REINVENT_WEB_JOBS_DIR"), repository_root / "jobs"
            ),
            app_timezone=os.getenv("REINVENT_WEB_TIMEZONE", "Asia/Tokyo"),
            cors_origins=origins,
            sqlite_busy_timeout_ms=int(
                os.getenv("REINVENT_WEB_SQLITE_BUSY_TIMEOUT_MS", "5000")
            ),
            worker_poll_interval_seconds=float(
                os.getenv("REINVENT_WEB_WORKER_POLL_SECONDS", "2.0")
            ),
            cancel_grace_seconds=float(
                os.getenv("REINVENT_WEB_CANCEL_GRACE_SECONDS", "5.0")
            ),
            reconciliation_grace_seconds=float(
                os.getenv("REINVENT_WEB_RECONCILIATION_GRACE_SECONDS", "30.0")
            ),
            parent_max_atoms=int(os.getenv("REINVENT_WEB_PARENT_MAX_ATOMS", "200")),
            worker_runner=os.getenv("REINVENT_WEB_WORKER_RUNNER", "fake"),
            reinvent_runtime_python=_resolve_path(
                os.getenv("REINVENT_RUNTIME_PYTHON"), Path(sys.executable)
            ),
            reinvent_source_root=_resolve_path(
                os.getenv("REINVENT_SOURCE_ROOT"), repository_root
            ),
        )


def get_settings() -> Settings:
    return Settings.from_env()
