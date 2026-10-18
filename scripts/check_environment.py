#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def check(label: str, python: Path, packages: list[str], command: list[str] | None = None) -> bool:
    if not python.is_file():
        print(f"[MISSING] {label}: {python}")
        return False
    script = "import importlib.metadata as m, json, sys; print(json.dumps({'python': sys.version.split()[0], 'packages': {p: m.version(p) for p in " + repr(packages) + "}}))"
    try:
        result = subprocess.run([str(python), "-c", script], check=True, capture_output=True, text=True, timeout=180)
        print(f"[OK] {label}: {result.stdout.strip()}")
        if command:
            subprocess.run([str(python), *command], check=True, capture_output=True, text=True, timeout=180)
        return True
    except (OSError, subprocess.SubprocessError) as error:
        print(f"[FAILED] {label}: {error}")
        return False


def main() -> int:
    reinvent_python = Path(os.environ.get("REINVENT_RUNTIME_PYTHON", ROOT / ".runtimes/reinvent/bin/python"))
    legacy_conda = Path.home() / "miniforge3/envs/r4/bin/python"
    if not reinvent_python.is_file() and legacy_conda.is_file():
        reinvent_python = legacy_conda
    runtimes = [
        ("REINVENT4 / RDKit / sklearn / boosters", reinvent_python, ["reinvent", "rdkit", "scikit-learn", "lightgbm", "xgboost", "catboost"], ["-m", "reinvent", "--help"]),
        ("Web backend", ROOT / "webapp/backend/.venv/bin/python", ["reinvent4-web-backend", "fastapi", "rdkit"], None),
        ("Chemprop v2", ROOT / ".runtimes/chemprop-v2/bin/python", ["chemprop", "torch", "rdkit"], ["-c", "import chemprop, torch, rdkit"]),
        ("QSARtuna", ROOT / ".runtimes/qsartuna/bin/python", ["qsartuna", "optuna", "chemprop", "peptides"], ["-c", "import os; os.environ['NUMBA_DISABLE_JIT']='1'; import qsartuna, peptides"]),
    ]
    results = [check(*item) for item in runtimes]
    print(json.dumps({"repository": str(ROOT), "all_required_checks_passed": all(results)}, indent=2))
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
