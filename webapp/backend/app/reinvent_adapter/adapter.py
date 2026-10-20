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
    raw_output: Path | None
    process_log: Path
    public_log: Path
    environment: dict[str, str]
    config_path: Path
    model_path: str
    output_model_path: Path | None = None


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
    resolved_model_path = str(conditions.get("resolved_model_path", model_spec.relative_path))
    repository_model_path = (settings.repository_root / resolved_model_path).resolve()
    source_model_path = (source_root / resolved_model_path).resolve()
    if repository_model_path.is_relative_to(settings.model_store_dir) and repository_model_path.is_file():
        model_path = repository_model_path
    elif source_model_path.is_relative_to(source_root) and source_model_path.is_file():
        model_path = source_model_path
    else:
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
    ]
    if model_spec.mode != "reinvent":
        lines.append(f"smiles_file = {_toml_string(input_file)}")
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
        # Keep the configured relative path in job metadata.  Custom TL/RL
        # models live under the repository model store, not under the
        # read-only REINVENT source tree, so ``relative_to(source_root)``
        # raises ValueError for those valid paths.
        model_path=resolved_model_path,
    )


def build_transfer_learning_command(
    settings: Settings, job: JobRecord, job_dir: Path
) -> ReinventCommand:
    conditions = job.run_conditions
    base_model_key = str(conditions.get("base_model_key", ""))
    if base_model_key not in {"reinvent", "mol2mol_scaffold_generic"}:
        raise AdapterConfigurationError("base_model_key is not supported for transfer learning")

    runtime_python = settings.reinvent_runtime_python.resolve()
    if not runtime_python.is_file() or not os.access(runtime_python, os.X_OK):
        raise AdapterConfigurationError("configured REINVENT Python is not executable")
    source_root = settings.reinvent_source_root.resolve()
    if not (source_root / "reinvent" / "__init__.py").is_file():
        raise AdapterConfigurationError("configured REINVENT source root is invalid")
    input_model = (source_root / MODEL_SPECS[base_model_key].relative_path).resolve()
    if not input_model.is_relative_to(source_root) or not input_model.is_file():
        raise AdapterConfigurationError("configured REINVENT base model is missing")

    output_relative = str(conditions.get("output_model_path", ""))
    output_model = (settings.repository_root / output_relative).resolve()
    if not output_model.is_relative_to(settings.model_store_dir) or output_model.suffix != ".model":
        raise AdapterConfigurationError("output model path is invalid")

    device = str(conditions.get("device", "mps"))
    if device not in {"cpu", "cuda", "mps"}:
        raise AdapterConfigurationError("device is not allowed")
    num_epochs = _bounded_int(
        conditions.get("num_epochs", 3), default=3, minimum=1, maximum=1000, name="num_epochs"
    )
    batch_size = _bounded_int(
        conditions.get("batch_size", 50), default=50, minimum=1, maximum=4096, name="batch_size"
    )
    learning_rate = _bounded_float(
        conditions.get("learning_rate", 1e-4),
        default=1e-4,
        minimum=1e-8,
        maximum=1.0,
        name="learning_rate",
    )
    lower = _bounded_float(
        conditions.get("tanimoto_lower", 0.7), default=0.7, minimum=0.0, maximum=1.0, name="tanimoto_lower"
    )
    upper = _bounded_float(
        conditions.get("tanimoto_upper", 1.0), default=1.0, minimum=0.0, maximum=1.0, name="tanimoto_upper"
    )
    if lower > upper:
        raise AdapterConfigurationError("tanimoto lower threshold exceeds upper threshold")

    train_file = job_dir / "input" / "train.smi"
    validation_file = job_dir / "input" / "validation.smi"
    public_log = job_dir / "logs" / "transfer-learning.log"
    process_log = job_dir / "logs" / "transfer-learning-console.log"
    json_config = job_dir / "config" / "transfer-learning.json"
    toml_path = job_dir / "config" / "transfer-learning.toml"
    tb_logdir = job_dir / "logs" / "tensorboard"
    progress_endpoint = f"{settings.internal_api_base_url}/api/training-jobs/{job.id}/progress"
    lines = [
        'run_type = "transfer_learning"',
        f"device = {_toml_string(device)}",
        f"tb_logdir = {_toml_string(tb_logdir)}",
        f"json_out_config = {_toml_string(json_config)}",
        "",
        "[parameters]",
        f"num_epochs = {num_epochs}",
        "save_every_n_epochs = 1",
        f"batch_size = {batch_size}",
        "num_refs = 0",
        "sample_batch_size = 100",
        f"input_model_file = {_toml_string(input_model)}",
        f"smiles_file = {_toml_string(train_file)}",
        f"validation_smiles_file = {_toml_string(validation_file)}",
        f"output_model_file = {_toml_string(output_model)}",
    ]
    if conditions.get("base_model") == "mol2mol":
        lines.extend(
            [
                'pairs.type = "tanimoto"',
                f"pairs.upper_threshold = {upper}",
                f"pairs.lower_threshold = {lower}",
                "pairs.min_cardinality = 1",
                "pairs.max_cardinality = 199",
            ]
        )
    lines.extend(
        [
            "",
            "[scheduler]",
            f"lr = {learning_rate}",
            "",
            "[responder]",
            f"endpoint = {_toml_string(progress_endpoint)}",
            "frequency = 1",
        ]
    )
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
        raw_output=None,
        process_log=process_log,
        public_log=public_log,
        environment=environment,
        config_path=toml_path,
        model_path=output_relative,
        output_model_path=output_model,
    )


def _rl_model_path(settings: Settings, value: str) -> Path:
    source_candidate = (settings.reinvent_source_root / value).resolve()
    repository_candidate = (settings.repository_root / value).resolve()
    if source_candidate.is_relative_to(settings.reinvent_source_root) and source_candidate.is_file():
        return source_candidate
    if repository_candidate.is_relative_to(settings.model_store_dir) and repository_candidate.is_file():
        return repository_candidate
    raise AdapterConfigurationError("configured RL model is missing or outside an allowed model directory")


def _transform_lines(transform: str, low: float, high: float, k: float) -> list[str]:
    if transform == "none":
        return []
    lines = [f'transform.type = {_toml_string(transform)}', f"transform.high = {high}", f"transform.low = {low}"]
    if transform == "double_sigmoid":
        lines.extend([f"transform.coef_div = {max(abs(high - low), 1.0)}", "transform.coef_si = 20.0", "transform.coef_se = 20.0"])
    else:
        lines.append(f"transform.k = {k}")
    return lines


def build_rl_command(settings: Settings, job: JobRecord, job_dir: Path) -> ReinventCommand:
    conditions = job.run_conditions
    runtime_python = settings.reinvent_runtime_python.resolve()
    source_root = settings.reinvent_source_root.resolve()
    if not runtime_python.is_file() or not os.access(runtime_python, os.X_OK):
        raise AdapterConfigurationError("configured REINVENT Python is not executable")
    if not (source_root / "reinvent" / "__init__.py").is_file():
        raise AdapterConfigurationError("configured REINVENT source root is invalid")
    prior_path = _rl_model_path(settings, str(conditions["prior_path"]))
    agent_path = _rl_model_path(settings, str(conditions["agent_path"]))
    output_model = (settings.repository_root / str(conditions["output_model_path"])).resolve()
    if not output_model.is_relative_to(settings.model_store_dir) or output_model.suffix != ".chkpt":
        raise AdapterConfigurationError("RL checkpoint path is invalid")
    summary_prefix = job_dir / "output" / "rl-summary"
    input_path = (settings.repository_root / str(conditions["input_smiles_path"])).resolve()
    config_path = job_dir / "config" / "reinforcement-learning.toml"
    json_config = job_dir / "config" / "reinforcement-learning.json"
    public_log = job_dir / "logs" / "reinforcement-learning.log"
    process_log = job_dir / "logs" / "reinforcement-learning-console.log"
    tb_logdir = job_dir / "logs" / "tensorboard"
    device = str(conditions["device"])
    if device not in {"cpu", "cuda", "mps"}:
        raise AdapterConfigurationError("RL device is not allowed")

    lines = [
        'run_type = "staged_learning"',
        f"device = {_toml_string(device)}",
        f"tb_logdir = {_toml_string(tb_logdir)}",
        f"json_out_config = {_toml_string(json_config)}",
        "",
        "[parameters]",
        f"prior_file = {_toml_string(prior_path)}",
        f"agent_file = {_toml_string(agent_path)}",
        f"summary_csv_prefix = {_toml_string(summary_prefix)}",
        "use_checkpoint = false",
        "purge_memories = false",
        f"batch_size = {int(conditions['batch_size'])}",
        "unique_sequences = true",
        "randomize_smiles = true",
    ]
    if conditions["model_type"] == "mol2mol":
        lines.extend([f"smiles_file = {_toml_string(input_path)}", 'sample_strategy = "multinomial"', "distance_threshold = 100"])
    lines.extend([
        "",
        "[learning_strategy]",
        'type = "dap"',
        f"sigma = {int(conditions['sigma'])}",
        f"rate = {float(conditions['learning_rate'])}",
    ])
    if conditions.get("diversity_filter"):
        lines.extend([
            "",
            "[diversity_filter]",
            'type = "IdenticalMurckoScaffold"',
            f"bucket_size = {int(conditions['diversity_bucket_size'])}",
            f"minscore = {float(conditions['diversity_min_score'])}",
            "minsimilarity = 0.4",
            "penalty_multiplier = 0.5",
        ])
    lines.extend([
        "",
        "[responder]",
        f"endpoint = {_toml_string(f'{settings.internal_api_base_url}/api/rl-jobs/{job.id}/progress')}",
        "frequency = 1",
        "",
        "[[stage]]",
        f"chkpt_file = {_toml_string(output_model)}",
        'termination = "simple"',
        f"max_score = {float(conditions['target_mean_score'])}",
        f"min_steps = {int(conditions['minimum_steps'])}",
        f"max_steps = {int(conditions['maximum_steps'])}",
        "",
        "[stage.scoring]",
        f"type = {_toml_string(conditions['aggregation'])}",
        "",
        "[[stage.scoring.component]]",
        "[stage.scoring.component.SmilesStandardizationEligibility]",
        "[[stage.scoring.component.SmilesStandardizationEligibility.endpoint]]",
        'name = "SMILES Eligibility"',
        'params.standardization_version = "1.0.0"',
    ])
    for component in conditions.get("qsar_components", []):
        model_path = (settings.repository_root / component["model_path"]).resolve()
        if not model_path.is_relative_to(settings.qsar_model_store_dir) or not model_path.is_file():
            raise AdapterConfigurationError("QSAR component model path is invalid")
        dataset_path = (settings.repository_root / component.get("dataset_path", "")).resolve()
        if not dataset_path.is_relative_to(settings.repository_root) or not dataset_path.is_file():
            raise AdapterConfigurationError("QSAR component training dataset is invalid")
        if component.get("algorithm") in {"chemprop_v2", "qsartuna"}:
            runtime = Path(component["external_runtime_python"]).resolve()
            predictor = Path(__file__).resolve().parents[2] / "worker" / "external_qsar_predictor.py"
            metadata_path = (settings.repository_root / component["evaluation_path"]).resolve()
            if not runtime.is_file() or not predictor.is_file() or not dataset_path.is_file() or not metadata_path.is_file():
                raise AdapterConfigurationError("External QSAR runtime or artifact is missing")
            arguments = " ".join([
                str(predictor), "--algorithm", component["algorithm"], "--model", str(model_path),
                "--training-dataset", str(dataset_path), "--ad-threshold",
                str(float(component["ad_threshold"]) if component["ad_behavior"] == "penalize" else 0.0),
                "--below-ad-score", "0.0", "--metadata", str(metadata_path),
            ])
            lines.extend([
                "", "[[stage.scoring.component]]", "[stage.scoring.component.ExternalProcess]",
                "[[stage.scoring.component.ExternalProcess.endpoint]]",
                f"name = {_toml_string('QSAR: ' + component['model_name'])}",
                f"weight = {float(component['weight'])}",
                f"params.executable = {_toml_string(runtime)}",
                f"params.args = {_toml_string(arguments)}",
            ])
            lines.extend(_transform_lines(component["transform"], float(component["low"]), float(component["high"]), float(component["k"])))
            if component["ad_behavior"] == "penalize":
                lines.extend([
                    "", "[[stage.scoring.component]]", "[stage.scoring.component.QsarApplicabilityDomain]",
                    "[[stage.scoring.component.QsarApplicabilityDomain.endpoint]]",
                    f"name = {_toml_string('QSAR AD: ' + component['model_name'])}",
                    'params.model_file = ""',
                    f"params.dataset_file = {_toml_string(dataset_path)}",
                    f"params.threshold = {float(component['ad_threshold'])}",
                ])
            continue
        lines.extend([
            "",
            "[[stage.scoring.component]]",
            "[stage.scoring.component.QsarJoblib]",
            "[[stage.scoring.component.QsarJoblib.endpoint]]",
            f"name = {_toml_string('QSAR: ' + component['model_name'])}",
            f"weight = {float(component['weight'])}",
            f"params.model_file = {_toml_string(model_path)}",
            f"params.dataset_file = {_toml_string((settings.repository_root / component['dataset_path']).resolve())}",
            f"params.task_type = {_toml_string(component['task_type'])}",
            f"params.ad_threshold = {float(component['ad_threshold']) if component['ad_behavior'] == 'penalize' else 0.0}",
            "params.below_ad_score = 0.0",
            f"params.class_index = {int(component['class_index'])}",
        ])
        lines.extend(_transform_lines(component["transform"], float(component["low"]), float(component["high"]), float(component["k"])))
        if component["ad_behavior"] == "penalize":
            lines.extend([
                "", "[[stage.scoring.component]]", "[stage.scoring.component.QsarApplicabilityDomain]",
                "[[stage.scoring.component.QsarApplicabilityDomain.endpoint]]",
                f"name = {_toml_string('QSAR AD: ' + component['model_name'])}",
                f"params.model_file = {_toml_string(model_path)}",
                f"params.dataset_file = {_toml_string(dataset_path)}",
                f"params.threshold = {float(component['ad_threshold'])}",
            ])
    for component in conditions.get("property_components", []):
        name = component["component"]
        lines.extend([
            "",
            "[[stage.scoring.component]]",
            f"[stage.scoring.component.{name}]",
            f"[[stage.scoring.component.{name}.endpoint]]",
            f"name = {_toml_string(name)}",
            f"weight = {float(component['weight'])}",
        ])
        if name == "TanimotoSimilarity":
            lines.extend([f"params.smiles = {json.dumps(component['reference_smiles'])}", "params.radius = 2", "params.use_counts = true", "params.use_features = true"])
        if name == "CustomAlerts":
            lines.append(f"params.smarts = {json.dumps(component['smarts'])}")
        lines.extend(_transform_lines(component["transform"], float(component["low"]), float(component["high"]), float(component["k"])))
    config_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    backend_root = Path(__file__).resolve().parents[2]
    return ReinventCommand(
        argv=[str(runtime_python), "-m", "reinvent", "-l", str(public_log), str(config_path)],
        raw_output=None,
        process_log=process_log,
        public_log=public_log,
        environment={"PYTHONPATH": os.pathsep.join((str(source_root), str(backend_root))), "PYTHONUNBUFFERED": "1"},
        config_path=config_path,
        model_path=str(conditions["output_model_path"]),
        output_model_path=output_model,
    )
