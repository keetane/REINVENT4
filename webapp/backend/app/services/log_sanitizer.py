from __future__ import annotations

import re
from pathlib import Path

from app.config import Settings


SECRET_ASSIGNMENT = re.compile(
    r"(?i)\b(api[_-]?key|access[_-]?token|auth[_-]?token|password|secret)\s*([=:])\s*([^\s,;]+)"
)


def sanitize_log_lines(lines: list[str], settings: Settings) -> list[str]:
    replacements = {
        str(Path.home().resolve()): "<home>",
        str(settings.repository_root.resolve()): "<repository>",
        str(settings.jobs_dir.resolve()): "<jobs>",
        str(settings.reinvent_source_root.resolve()): "<reinvent-source>",
        str(settings.reinvent_runtime_python.resolve()): "<reinvent-python>",
    }
    ordered = sorted(replacements.items(), key=lambda item: len(item[0]), reverse=True)
    sanitized: list[str] = []
    for line in lines:
        current = SECRET_ASSIGNMENT.sub(lambda match: f"{match.group(1)}{match.group(2)}<redacted>", line)
        for raw, replacement in ordered:
            if raw:
                current = current.replace(raw, replacement)
        sanitized.append(current)
    return sanitized
