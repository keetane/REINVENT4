from __future__ import annotations

import re
import unicodedata
from datetime import datetime
from pathlib import Path
from uuid import UUID
from zoneinfo import ZoneInfo


SEPARATOR_RE = re.compile(r"[_.-]+")

MODEL_TYPE_LABELS = {
    "reinvent": "Reinvent",
    "libinvent": "LibInvent",
    "linkinvent": "LinkInvent",
    "mol2mol_high": "Mol2MolHigh",
    "mol2mol_medium": "Mol2MolMedium",
    "mol2mol_mmp": "Mol2MolMMP",
    "mol2mol_scaffold": "Mol2MolScaffold",
    "mol2mol_scaffold_generic": "Mol2MolScaffoldGeneric",
    "mol2mol_similarity": "Mol2MolSimilarity",
    "chemprop_v2": "ChempropV2",
    "qsartuna": "QSARtuna",
    "lightgbm": "LightGBM",
    "xgboost": "XGBoost",
    "catboost": "CatBoost",
    "mlp": "MLP",
    "random_forest": "RandomForest",
    "svr": "SVR",
    "elastic_net": "ElasticNet",
}


def slugify_comment(comment: str, *, ascii_only: bool = False, max_length: int = 50) -> str:
    normalized = unicodedata.normalize("NFKC", comment).strip()
    output: list[str] = []
    for character in normalized:
        category = unicodedata.category(character)
        if character.isspace():
            output.append("_")
        elif character in {"_", "-", "."}:
            output.append(character)
        elif category.startswith("C") or character in {"/", "\\"}:
            output.append("_")
        elif character.isalnum() and (not ascii_only or character.isascii()):
            output.append(character)
    slug = SEPARATOR_RE.sub(lambda match: match.group(0)[0], "".join(output))
    slug = slug.lstrip(".").strip("_.-")[:max_length].rstrip("_.-")
    return slug if slug not in {"", ".", ".."} else "result"


def result_filename(
    *, comment_slug: str, created_at: datetime, timezone_name: str,
    job_id: UUID | None = None, tab_type: str | None = None,
    model_type: str | None = None,
) -> str:
    local_time = created_at.astimezone(ZoneInfo(timezone_name))
    suffix = f"_{str(job_id)[:8]}" if job_id is not None else ""
    if tab_type and model_type:
        return f"{local_time:%Y%m%d}_{tab_type}_{model_type}_{comment_slug}{suffix}.csv"
    return f"{local_time:%Y%m%d%H%M}_{comment_slug}{suffix}.csv"


def model_type_label(model_key: str) -> str:
    return MODEL_TYPE_LABELS.get(model_key, slugify_comment(model_key, ascii_only=True))


def next_artifact_filename(
    directory: Path,
    *,
    tab_type: str,
    model_key: str,
    comment: str,
    date: datetime,
    extension: str,
) -> str:
    comment_slug = slugify_comment(comment)
    base = f"{date:%Y%m%d}_{tab_type}_{model_type_label(model_key)}_{comment_slug}"
    suffix = "" if not extension else extension if extension.startswith(".") else f".{extension}"
    candidate = f"{base}{suffix}"
    serial = 2
    def occupied(value: str) -> bool:
        if suffix:
            return (directory / value).exists()
        return any(path.stem == value for path in directory.iterdir())

    while occupied(candidate):
        candidate = f"{base}_{serial}{suffix}"
        serial += 1
    return candidate


def confined_path(root: Path, relative_path: str) -> Path:
    root_resolved = root.resolve()
    candidate = (root_resolved / relative_path).resolve()
    if not candidate.is_relative_to(root_resolved):
        raise ValueError("path escapes configured repository root")
    return candidate
