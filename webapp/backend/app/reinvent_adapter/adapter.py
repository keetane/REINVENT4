from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from app.config import Settings
from app.db.jobs import JobRecord
from app.reinvent_adapter.models import MODEL_SPECS


class AdapterConfigurationError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class ReinventCommand:
    argv: list[str]
    raw_output: Path
    process_log: Path
    public_log: Path
    environment: dict[str, str]
    config_path: Path
    model_path: str


def _toml_string(value: str | Path) -> str:
    return json.dumps(str(value), ensure_ascii=False)


def _bounded_int(value, *, default: int, minimum: int, maximum: int, name: str) -> int:
    if isinstance(value, bool):
        raise AdapterConfigurationError(f"{name} must be an integer")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise AdapterConfigurationError(f"{name} must be an integer") from error
    if not minimum <= parsed <= maximum:
        raise AdapterConfigurationError(
            f"{name} must be between {minimum} and {maximum}"
        )
    return parsed


def _boolean(value, *, name: str) -> bool:
    if not isinstance(value, bool):
        raise AdapterConfigurationError(f"{name} must be a boolean")
    return value


def _bounded_float(
    value, *, default: float, minimum: float, maximum: float, name: str
) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as error:
        raise AdapterConfigurationError(f"{name} must be numeric") from error
    if not minimum <= parsed <= maximum:
        raise AdapterConfigurationError(
            f"{name} must be between {minimum} and {maximum}"
        )
    return parsed


def build_reinvent_command(
    settings: Settings, job: JobRecord, job_dir: Path
) -> ReinventCommand:
    conditions = job.run_conditions
    model_key = str(conditions.get("model_key", "mol2mol_high"))
    try:
        model_spec = MODEL_SPECS[model_key]
    except KeyError as error:
        raise AdapterConfigurationError("model_key is not in the server allowlist") from error

    runtime_python = settings.reinvent_runtime_python.resolve()
    if not runtime_python.is_file() or not os.access(runtime_python, os.X_OK):
        raise AdapterConfigurationError("configured REINVENT Python is not executable")
    source_root = settings.reinvent_source_root.resolve()
    if not (source_root / "reinvent" / "__init__.py").is_file():
        raise AdapterConfigurationError("configured REINVENT source root is invalid")
    model_path = (source_root / model_spec.relative_path).resolve()
    if not model_path.is_relative_to(source_root) or not model_path.is_file():
        raise AdapterConfigurationError("configured REINVENT model is missing")

    num_smiles = _bounded_int(
        conditions.get("num_smiles", 100),
        default=100,
        minimum=1,
        maximum=100_000,
        name="num_smiles",
    )
    device = str(conditions.get("device", "cpu"))
    if device not in {"cpu", "cuda", "mps"}:
        raise AdapterConfigurationError("device is not allowed")
    strategy = str(conditions.get("sample_strategy", "multinomial"))
    if strategy not in {"beamsearch", "multinomial"}:
        raise AdapterConfigurationError("sample_strategy is not allowed")
    temperature = _bounded_float(
        conditions.get("temperature", 1.0),
        default=1.0,
        minimum=0.0,
        maximum=1.0,
        name="temperature",
    )
    unique_molecules = _boolean(
        conditions.get("unique_molecules", True), name="unique_molecules"
    )
    randomize_smiles = _boolean(
        conditions.get("randomize_smiles", True), name="randomize_smiles"
    )

    input_file = job_dir / "input" / "parent.smi"
    raw_output = job_dir / "output" / "reinvent-raw.csv"
    public_log = job_dir / "logs" / "reinvent.log"
    process_log = job_dir / "logs" / "reinvent-console.log"
    json_config = job_dir / "config" / "sampling.json"
    toml_path = job_dir / "config" / "sampling.toml"

    lines = [
        'run_type = "sampling"',
        f"device = {_toml_string(device)}",
        f"json_out_config = {_toml_string(json_config)}",
        "",
        "[parameters]",
        f"model_file = {_toml_string(model_path)}",
        f"output_file = {_toml_string(raw_output)}",
        f"num_smiles = {num_smiles}",
        f"unique_molecules = {str(unique_molecules).lower()}",
        f"randomize_smiles = {str(randomize_smiles).lower()}",
        f"smiles_file = {_toml_string(input_file)}",
    ]
    if model_spec.mode == "mol2mol":
        lines.append(f"sample_strategy = {_toml_string(strategy)}")
        if strategy == "multinomial":
            lines.append(f"temperature = {temperature}")
    toml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    environment = {"PYTHONPATH": str(source_root), "PYTHONUNBUFFERED": "1"}
    return ReinventCommand(
        argv=[
            str(runtime_python),
            "-m",
            "reinvent",
            "-l",
            str(public_log),
            str(toml_path),
        ],
        raw_output=raw_output,
        process_log=process_log,
        public_log=public_log,
        environment=environment,
        config_path=toml_path,
        model_path=str(model_path.relative_to(source_root)),
    )
