from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any
from uuid import UUID

from app.config import Settings
from app.db.connection import connect
from app.db.jobs import utc_now


@dataclass(frozen=True, slots=True)
class TrainedModelRecord:
    id: UUID
    name: str
    model_type: str
    created_at: str
    training_data_count: int
    validation_data_count: int
    training_conditions: dict[str, Any]
    file_path: str
    config_path: str
    log_path: str
    source_job_id: UUID
    base_model_key: str


def _record(row) -> TrainedModelRecord:
    return TrainedModelRecord(
        id=UUID(row["id"]),
        name=row["name"],
        model_type=row["model_type"],
        created_at=row["created_at"],
        training_data_count=row["training_data_count"],
        validation_data_count=row["validation_data_count"],
        training_conditions=json.loads(row["training_conditions_json"]),
        file_path=row["file_path"],
        config_path=row["config_path"],
        log_path=row["log_path"],
        source_job_id=UUID(row["source_job_id"]),
        base_model_key=row["base_model_key"],
    )


class ModelRepository:
    def __init__(self, settings: Settings):
        self.settings = settings

    def register(
        self,
        *,
        model_id: UUID,
        name: str,
        model_type: str,
        training_data_count: int,
        validation_data_count: int,
        training_conditions: dict[str, Any],
        file_path: str,
        config_path: str,
        log_path: str,
        source_job_id: UUID,
        base_model_key: str,
    ) -> TrainedModelRecord:
        with connect(self.settings) as connection:
            connection.execute(
                """
                INSERT INTO trained_models (
                    id, name, model_type, created_at, training_data_count,
                    validation_data_count, training_conditions_json, file_path,
                    config_path, log_path, source_job_id, base_model_key
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    str(model_id),
                    name,
                    model_type,
                    utc_now(),
                    training_data_count,
                    validation_data_count,
                    json.dumps(training_conditions, ensure_ascii=False, sort_keys=True),
                    file_path,
                    config_path,
                    log_path,
                    str(source_job_id),
                    base_model_key,
                ),
            )
        return self.get(model_id)

    def get(self, model_id: UUID) -> TrainedModelRecord:
        with connect(self.settings) as connection:
            row = connection.execute(
                "SELECT * FROM trained_models WHERE id = ?", (str(model_id),)
            ).fetchone()
        if row is None:
            raise KeyError(str(model_id))
        return _record(row)

    def list(self) -> list[TrainedModelRecord]:
        with connect(self.settings) as connection:
            rows = connection.execute(
                "SELECT * FROM trained_models ORDER BY created_at DESC"
            ).fetchall()
        return [_record(row) for row in rows]

    def is_in_use(self, model_id: UUID) -> bool:
        with connect(self.settings) as connection:
            row = connection.execute(
                """
                SELECT 1 FROM jobs
                WHERE status IN ('queued', 'running', 'cancelling')
                    AND (model_id = ? OR run_conditions_json LIKE ?)
                LIMIT 1
                """,
                (str(model_id), f'%"{model_id}"%'),
            ).fetchone()
        return row is not None

    def delete(self, model_id: UUID) -> None:
        with connect(self.settings) as connection:
            deleted = connection.execute(
                "DELETE FROM trained_models WHERE id = ?", (str(model_id),)
            ).rowcount
        if deleted != 1:
            raise KeyError(str(model_id))
