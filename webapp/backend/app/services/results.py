from __future__ import annotations

import csv
import os
from pathlib import Path


MAX_RESULT_BYTES = 100 * 1024 * 1024


class ResultValidationError(ValueError):
    pass


def validate_and_finalize_csv(source: Path, destination: Path, job_dir: Path) -> int:
    job_root = job_dir.resolve()
    source_resolved = source.resolve()
    destination_resolved = destination.resolve()
    if not source_resolved.is_relative_to(job_root):
        raise ResultValidationError("result source escapes job directory")
    if not destination_resolved.is_relative_to(job_root):
        raise ResultValidationError("result destination escapes job directory")
    if source.is_symlink() or not source.is_file():
        raise ResultValidationError("expected result CSV is missing or is not a regular file")
    size = source.stat().st_size
    if size <= 0:
        raise ResultValidationError("result CSV is empty")
    if size > MAX_RESULT_BYTES:
        raise ResultValidationError("result CSV exceeds the configured size limit")
    try:
        with source.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "SMILES" not in reader.fieldnames:
                raise ResultValidationError("result CSV must contain a SMILES column")
            rows = sum(1 for _ in reader)
    except UnicodeDecodeError as error:
        raise ResultValidationError("result CSV is not valid UTF-8") from error
    if rows == 0:
        raise ResultValidationError("result CSV contains no data rows")
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, destination)
    return rows
