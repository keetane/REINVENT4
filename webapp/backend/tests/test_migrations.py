from app.db.connection import connect
from app.db.migrations import apply_migrations, current_versions


def test_migration_is_idempotent_and_enables_required_pragmas(settings):
    assert apply_migrations(settings) == ["0001_create_jobs", "0002_transfer_learning", "0003_qsar_models", "0004_expand_qsar_algorithms"]
    assert apply_migrations(settings) == []
    assert current_versions(settings) == ["0001_create_jobs", "0002_transfer_learning", "0003_qsar_models", "0004_expand_qsar_algorithms"]

    with connect(settings) as connection:
        assert connection.execute("PRAGMA foreign_keys").fetchone()[0] == 1
        assert connection.execute("PRAGMA journal_mode").fetchone()[0] == "wal"
        columns = {
            row["name"] for row in connection.execute("PRAGMA table_info(jobs)")
        }

    assert {"id", "status", "job_dir", "run_conditions_json", "job_type", "model_id"} <= columns
    with connect(settings) as connection:
        sql = connection.execute(
            "SELECT sql FROM sqlite_master WHERE type='table' AND name='jobs'"
        ).fetchone()[0]
        assert "qsar_training" in sql
        assert connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='qsar_models'"
        ).fetchone() is not None
