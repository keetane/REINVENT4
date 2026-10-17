from __future__ import annotations

import json
from pathlib import Path
from uuid import UUID

from app.config import Settings


def relative_to_repository(settings: Settings, path: Path) -> str:
    return str(path.resolve().relative_to(settings.repository_root.resolve()))


def create_job_workspace(
    settings: Settings,
    *,
    job_id: UUID,
    parent_smiles: str,
    comment: str,
    run_conditions: dict,
) -> Path:
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    if not job_dir.is_relative_to(settings.jobs_dir.resolve()):
        raise ValueError("job directory escaped configured jobs root")
    for name in ("input", "config", "logs", "output"):
        (job_dir / name).mkdir(parents=True, exist_ok=False)
    (job_dir / "input" / "parent.smi").write_text(parent_smiles + "\n", encoding="utf-8")
    snapshot = {
        "comment": comment,
        "parent_smiles": parent_smiles,
        "run_conditions": run_conditions,
    }
    (job_dir / "input" / "request.json").write_text(
        json.dumps(snapshot, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return job_dir

