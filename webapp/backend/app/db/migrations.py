from __future__ import annotations

from pathlib import Path

from app.config import Settings
from app.db.connection import connect


MIGRATIONS_DIR = Path(__file__).with_name("migrations")


def apply_migrations(settings: Settings) -> list[str]:
    applied_now: list[str] = []
    with connect(settings) as connection:
        connection.execute(
            """
            CREATE TABLE IF NOT EXISTS schema_migrations (
                version TEXT PRIMARY KEY,
                applied_at TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        applied = {
            row["version"]
            for row in connection.execute("SELECT version FROM schema_migrations")
        }
        for migration_path in sorted(MIGRATIONS_DIR.glob("*.sql")):
            version = migration_path.stem
            if version in applied:
                continue
            escaped_version = version.replace("'", "''")
            script = migration_path.read_text(encoding="utf-8")
            try:
                connection.executescript(
                    "BEGIN IMMEDIATE;\n"
                    f"{script}\n"
                    "INSERT INTO schema_migrations(version) "
                    f"VALUES ('{escaped_version}');\n"
                    "COMMIT;"
                )
            except Exception:
                connection.execute("ROLLBACK")
                raise
            applied_now.append(version)
    return applied_now


def current_versions(settings: Settings) -> list[str]:
    with connect(settings) as connection:
        row = connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='schema_migrations'"
        ).fetchone()
        if row is None:
            return []
        return [
            item["version"]
            for item in connection.execute(
                "SELECT version FROM schema_migrations ORDER BY version"
            )
        ]

