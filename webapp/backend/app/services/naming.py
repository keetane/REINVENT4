from __future__ import annotations

import re
import unicodedata
from datetime import datetime
from pathlib import Path
from uuid import UUID
from zoneinfo import ZoneInfo


SEPARATOR_RE = re.compile(r"[_.-]+")


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
    *, comment_slug: str, created_at: datetime, timezone_name: str, job_id: UUID | None = None
) -> str:
    local_time = created_at.astimezone(ZoneInfo(timezone_name))
    suffix = f"_{str(job_id)[:8]}" if job_id is not None else ""
    return f"{local_time:%Y%m%d%H%M}_{comment_slug}{suffix}.csv"


def confined_path(root: Path, relative_path: str) -> Path:
    root_resolved = root.resolve()
    candidate = (root_resolved / relative_path).resolve()
    if not candidate.is_relative_to(root_resolved):
        raise ValueError("path escapes configured repository root")
    return candidate

