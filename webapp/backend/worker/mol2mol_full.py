from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path


SOURCE_MODEL_COLUMN = "GenerationModel"
ALLOWED_MODEL_KEYS = {
    "mol2mol_high",
    "mol2mol_medium",
    "mol2mol_mmp",
    "mol2mol_scaffold",
    "mol2mol_scaffold_generic",
}


def _confined(path_value: str, root: Path) -> Path:
    path = Path(path_value).resolve()
    if not path.is_relative_to(root):
        raise ValueError("Mol2Mol Full manifest contains a path outside the job directory")
    return path


def combine_results(runs: list[dict[str, str]], destination: Path, job_dir: Path) -> int:
    """Combine model outputs while keeping GenerationModel as the final column."""
    fieldnames: list[str] = []
    sources: list[tuple[str, Path]] = []
    for run in runs:
        model_key = str(run.get("model_key", ""))
        if model_key not in ALLOWED_MODEL_KEYS:
            raise ValueError("Mol2Mol Full manifest contains an unsupported model")
        source = _confined(str(run.get("raw_output", "")), job_dir)
        if source.is_symlink() or not source.is_file():
            raise ValueError(f"Mol2Mol output is missing: {model_key}")
        with source.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None or "SMILES" not in reader.fieldnames:
                raise ValueError(f"Mol2Mol output has no SMILES column: {model_key}")
            for name in reader.fieldnames:
                if name != SOURCE_MODEL_COLUMN and name not in fieldnames:
                    fieldnames.append(name)
        sources.append((model_key, source))

    destination = _confined(str(destination), job_dir)
    destination.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=[*fieldnames, SOURCE_MODEL_COLUMN], extrasaction="ignore"
        )
        writer.writeheader()
        for model_key, source in sources:
            with source.open(newline="", encoding="utf-8") as source_handle:
                for row in csv.DictReader(source_handle):
                    writer.writerow({**row, SOURCE_MODEL_COLUMN: model_key})
                    row_count += 1
    if row_count == 0:
        raise ValueError("Mol2Mol Full produced no molecules")
    return row_count


def run_manifest(manifest_path: Path) -> int:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    job_dir = Path(str(manifest["job_dir"])).resolve()
    if not manifest_path.resolve().is_relative_to(job_dir):
        raise ValueError("Mol2Mol Full manifest is outside the job directory")
    runtime_python = Path(str(manifest["runtime_python"])).resolve()
    source_root = Path(str(manifest["source_root"])).resolve()
    if not runtime_python.is_file() or not os.access(runtime_python, os.X_OK):
        raise ValueError("configured REINVENT Python is not executable")
    if not (source_root / "reinvent" / "__init__.py").is_file():
        raise ValueError("configured REINVENT source root is invalid")
    runs = manifest.get("runs")
    if (
        not isinstance(runs, list)
        or len(runs) != len(ALLOWED_MODEL_KEYS)
        or {str(run.get("model_key", "")) for run in runs if isinstance(run, dict)}
        != ALLOWED_MODEL_KEYS
    ):
        raise ValueError("Mol2Mol Full manifest must contain all five models")

    environment = dict(os.environ)
    environment.update({"PYTHONPATH": str(source_root), "PYTHONUNBUFFERED": "1"})
    for index, run in enumerate(runs, start=1):
        if not isinstance(run, dict):
            raise ValueError("Mol2Mol Full manifest contains an invalid run")
        model_key = str(run.get("model_key", ""))
        if model_key not in ALLOWED_MODEL_KEYS:
            raise ValueError("Mol2Mol Full manifest contains an unsupported model")
        config_path = _confined(str(run.get("config_path", "")), job_dir)
        public_log = _confined(str(run.get("public_log", "")), job_dir)
        print(f"[{index}/5] Sampling with {model_key}", flush=True)
        result = subprocess.run(
            [
                str(runtime_python),
                "-m",
                "reinvent",
                "-l",
                str(public_log),
                str(config_path),
            ],
            cwd=job_dir,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=sys.stdout,
            stderr=subprocess.STDOUT,
            check=False,
        )
        if result.returncode != 0:
            print(f"{model_key} failed with exit code {result.returncode}", flush=True)
            return result.returncode

    rows = combine_results(
        runs,
        _confined(str(manifest["combined_output"]), job_dir),
        job_dir,
    )
    print(f"Mol2Mol Full combined {rows} molecules", flush=True)
    return 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    try:
        raise SystemExit(run_manifest(args.manifest.resolve()))
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        print(f"Mol2Mol Full failed: {error}", flush=True)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
