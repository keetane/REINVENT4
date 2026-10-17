from __future__ import annotations

import argparse
import logging

from app.config import Settings, get_settings
from app.db.migrations import apply_migrations
from worker.engine import Worker, reconcile_active_jobs


LOGGER = logging.getLogger("reinvent.worker")


def initialize_worker(settings: Settings) -> None:
    settings.jobs_dir.mkdir(parents=True, exist_ok=True)
    apply_migrations(settings)
    reconcile_active_jobs(settings)


def run(settings: Settings, *, once: bool, poll_interval: float) -> None:
    initialize_worker(settings)
    worker = Worker(settings)
    LOGGER.info("worker ready: %s", worker.worker_id)
    if once:
        worker.run_once()
        return
    if poll_interval != settings.worker_poll_interval_seconds:
        settings = Settings(
            repository_root=settings.repository_root,
            database_path=settings.database_path,
            jobs_dir=settings.jobs_dir,
            app_timezone=settings.app_timezone,
            cors_origins=settings.cors_origins,
            sqlite_busy_timeout_ms=settings.sqlite_busy_timeout_ms,
            worker_poll_interval_seconds=poll_interval,
            cancel_grace_seconds=settings.cancel_grace_seconds,
            reconciliation_grace_seconds=settings.reconciliation_grace_seconds,
            parent_max_atoms=settings.parent_max_atoms,
            worker_runner=settings.worker_runner,
            reinvent_runtime_python=settings.reinvent_runtime_python,
            reinvent_source_root=settings.reinvent_source_root,
        )
        worker = Worker(settings, worker_id=worker.worker_id)
    worker.run_forever()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the independent REINVENT4 worker")
    parser.add_argument("--once", action="store_true", help="initialize and exit")
    parser.add_argument("--poll-interval", type=float, default=2.0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    run(get_settings(), once=args.once, poll_interval=args.poll_interval)


if __name__ == "__main__":
    main()
