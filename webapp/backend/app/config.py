from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]


def _resolve_path(value: str | None, default: Path) -> Path:
    return Path(value).expanduser().resolve() if value else default.resolve()


def _executable_path(value: str | None, default: Path) -> Path:
    """Return an absolute executable path without dereferencing virtualenv links.

    A virtual environment's ``bin/python`` is commonly a symlink.  Calling
    ``resolve()`` here changes the executable's location and therefore its
    ``sys.prefix``, which silently escapes the selected environment.
    """
    return Path(os.path.abspath(Path(value).expanduser() if value else default))


@dataclass(frozen=True, slots=True)
class Settings:
    repository_root: Path
    database_path: Path
    jobs_dir: Path
    models_dir: Path | None = None
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
    internal_api_base_url: str = "http://127.0.0.1:8000"
    chemprop_runtime_python: Path | None = None
    qsartuna_runtime_python: Path | None = None

    @property
    def model_store_dir(self) -> Path:
        return (self.models_dir or (self.repository_root / "models")).resolve()

    @property
    def qsar_model_store_dir(self) -> Path:
        return (self.model_store_dir / "qsar").resolve()

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
            models_dir=_resolve_path(
                os.getenv("REINVENT_WEB_MODELS_DIR"), repository_root / "models"
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
            reinvent_runtime_python=_executable_path(
                os.getenv("REINVENT_RUNTIME_PYTHON"), Path(sys.executable)
            ),
            reinvent_source_root=_resolve_path(
                os.getenv("REINVENT_SOURCE_ROOT"), repository_root
            ),
            internal_api_base_url=os.getenv(
                "REINVENT_WEB_INTERNAL_API_BASE_URL", "http://127.0.0.1:8000"
            ).rstrip("/"),
            chemprop_runtime_python=_executable_path(
                os.getenv("CHEMPROP_V2_PYTHON"), repository_root / ".runtimes" / "chemprop-v2" / "bin" / "python"
            ),
            qsartuna_runtime_python=_executable_path(
                os.getenv("QSARTUNA_PYTHON"), repository_root / ".runtimes" / "qsartuna" / "bin" / "python"
            ),
        )


def get_settings() -> Settings:
    return Settings.from_env()
