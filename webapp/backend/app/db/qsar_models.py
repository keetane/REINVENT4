from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any
from uuid import UUID

from app.config import Settings
from app.db.connection import connect
from app.db.jobs import utc_now


@dataclass(frozen=True, slots=True)
class QsarModelRecord:
    id: UUID
    name: str
    algorithm: str
    task_type: str
    target_name: str
    target_unit: str
    created_at: str
    training_data_count: int
    test_data_count: int
    metrics: dict[str, Any]
    metadata: dict[str, Any]
    file_path: str
    dataset_path: str
    predictions_path: str
    evaluation_path: str
    source_job_id: UUID
    ad_threshold: float


def _record(row) -> QsarModelRecord:
    return QsarModelRecord(
        id=UUID(row["id"]),
        name=row["name"],
        algorithm=row["algorithm"],
        task_type=row["task_type"],
        target_name=row["target_name"],
        target_unit=row["target_unit"],
        created_at=row["created_at"],
        training_data_count=row["training_data_count"],
        test_data_count=row["test_data_count"],
        metrics=json.loads(row["metrics_json"]),
        metadata=json.loads(row["metadata_json"]),
        file_path=row["file_path"],
        dataset_path=row["dataset_path"],
        predictions_path=row["predictions_path"],
        evaluation_path=row["evaluation_path"],
        source_job_id=UUID(row["source_job_id"]),
        ad_threshold=float(row["ad_threshold"]),
    )


class QsarModelRepository:
    def __init__(self, settings: Settings):
        self.settings = settings

    def register(
        self,
        *,
        model_id: UUID,
        name: str,
        algorithm: str,
        task_type: str,
        target_name: str,
        target_unit: str,
        training_data_count: int,
        test_data_count: int,
        metrics: dict[str, Any],
        metadata: dict[str, Any],
        file_path: str,
        dataset_path: str,
        predictions_path: str,
        evaluation_path: str,
        source_job_id: UUID,
        ad_threshold: float,
    ) -> QsarModelRecord:
        with connect(self.settings) as connection:
            connection.execute(
                """
                INSERT INTO qsar_models (
                    id, name, algorithm, task_type, target_name, target_unit,
                    created_at, training_data_count, test_data_count,
                    metrics_json, metadata_json, file_path, dataset_path,
                    predictions_path, evaluation_path, source_job_id, ad_threshold
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    str(model_id),
                    name,
                    algorithm,
                    task_type,
                    target_name,
                    target_unit,
                    utc_now(),
                    training_data_count,
                    test_data_count,
                    json.dumps(metrics, ensure_ascii=False, sort_keys=True),
                    json.dumps(metadata, ensure_ascii=False, sort_keys=True),
                    file_path,
                    dataset_path,
                    predictions_path,
                    evaluation_path,
                    str(source_job_id),
                    ad_threshold,
                ),
            )
        return self.get(model_id)

    def get(self, model_id: UUID) -> QsarModelRecord:
        with connect(self.settings) as connection:
            row = connection.execute(
                "SELECT * FROM qsar_models WHERE id = ?", (str(model_id),)
            ).fetchone()
        if row is None:
            raise KeyError(str(model_id))
        return _record(row)

    def list(self) -> list[QsarModelRecord]:
        with connect(self.settings) as connection:
            rows = connection.execute(
                "SELECT * FROM qsar_models ORDER BY created_at DESC"
            ).fetchall()
        return [_record(row) for row in rows]

    def is_in_use(self, model_id: UUID) -> bool:
        with connect(self.settings) as connection:
            row = connection.execute(
                """
                SELECT 1 FROM jobs
                WHERE job_type = 'reinforcement_learning'
                    AND status IN ('queued', 'running', 'cancelling')
                    AND run_conditions_json LIKE ?
                LIMIT 1
                """,
                (f'%"{model_id}"%',),
            ).fetchone()
        return row is not None

    def delete(self, model_id: UUID) -> None:
        with connect(self.settings) as connection:
            deleted = connection.execute(
                "DELETE FROM qsar_models WHERE id = ?", (str(model_id),)
            ).rowcount
        if deleted != 1:
            raise KeyError(str(model_id))
