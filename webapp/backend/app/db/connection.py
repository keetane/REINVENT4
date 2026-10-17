from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from collections.abc import Iterator

from app.config import Settings


@contextmanager
def connect(settings: Settings) -> Iterator[sqlite3.Connection]:
    settings.database_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(
        settings.database_path,
        timeout=settings.sqlite_busy_timeout_ms / 1_000,
        isolation_level=None,
    )
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA journal_mode = WAL")
    connection.execute(f"PRAGMA busy_timeout = {settings.sqlite_busy_timeout_ms}")
    try:
        yield connection
    finally:
        connection.close()

