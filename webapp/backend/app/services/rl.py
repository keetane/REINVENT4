from __future__ import annotations

import json
import math
import os
import subprocess
import csv
from functools import lru_cache
from datetime import datetime
from pathlib import Path
from uuid import UUID, uuid4
from zoneinfo import ZoneInfo

from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Contrib.SA_Score import sascorer

from app.config import Settings
from app.db.jobs import JobRecord, JobRepository
from app.db.models import ModelRepository
from app.db.qsar_models import QsarModelRepository
from app.errors import AppError
from app.reinvent_adapter.models import MODEL_SPECS
from app.schemas.rl import QsarPredictionRequest, RlJobCreateRequest, RlPreviewRequest
from app.services.naming import next_artifact_filename
from app.services.workspaces import relative_to_repository
from app.services.smiles_standardization import standardization_counts, standardization_metadata, standardize_smiles
from app.services.training import checkpoint_compatibility


@lru_cache(maxsize=4)
def _torch_devices(runtime_python: str) -> dict[str, bool]:
    code = (
        "import json,torch;"
        "print(json.dumps({'cpu':True,'cuda':torch.cuda.is_available(),"
        "'mps':bool(hasattr(torch.backends,'mps') and torch.backends.mps.is_available())}))"
    )
    try:
        completed = subprocess.run(
            [runtime_python, "-c", code], check=True, capture_output=True, text=True, timeout=20
        )
        return json.loads(completed.stdout)
    except (OSError, subprocess.SubprocessError, json.JSONDecodeError):
        return {"cpu": True, "cuda": False, "mps": False}


def available_torch_devices(settings: Settings) -> dict[str, bool]:
    return _torch_devices(str(settings.reinvent_runtime_python))


def _generation_model(settings: Settings, model_id: str):
    if model_id.startswith("standard:"):
        key = model_id.removeprefix("standard:")
        if key not in {"reinvent", "mol2mol_scaffold_generic"}:
            raise AppError(status_code=422, code="unsupported_rl_model", message="RLはREINVENTまたはMol2Molモデルに対応します。")
        spec = MODEL_SPECS[key]
        return key, spec.mode, spec.relative_path, spec.relative_path
    try:
        model = ModelRepository(settings).get(UUID(model_id))
    except (ValueError, KeyError) as error:
        raise AppError(status_code=404, code="generation_model_not_found", message="生成モデルが見つかりません。") from error
    prior_key = model.base_model_key
    if prior_key not in MODEL_SPECS or model.model_type not in {"reinvent", "mol2mol"}:
        raise AppError(status_code=422, code="incompatible_rl_model", message="この生成モデルはRLに対応していません。")
    return prior_key, model.model_type, MODEL_SPECS[prior_key].relative_path, model.file_path


def _next_rl_model_name(settings: Settings, model_key: str, comment: str) -> str:
    return next_artifact_filename(
        settings.model_store_dir, tab_type="RL", model_key=model_key,
        comment=comment, date=datetime.now(ZoneInfo(settings.app_timezone)), extension=".chkpt"
    )


def _input_smiles(settings: Settings, request: RlJobCreateRequest, model_type: str, agent_path: str) -> tuple[list[str], list, tuple[tuple[str, str], ...]]:
    values: list[str] = []
    results = []
    for line in request.input_smiles_content.splitlines():
        raw = line.strip().split(",", 1)[0].strip()
        if not raw:
            continue
        result = standardize_smiles(raw)
        results.append(result)
        if result.accepted and result.standardized_smiles and result.standardized_smiles not in values:
            values.append(result.standardized_smiles)
    source_candidate = (settings.reinvent_source_root / agent_path).resolve()
    repository_candidate = (settings.repository_root / agent_path).resolve()
    model_path = source_candidate if source_candidate.is_file() else repository_candidate
    values, compatibility_rejections = checkpoint_compatibility(settings, model_path, values)
    if model_type == "mol2mol" and not values:
        raise AppError(status_code=422, code="rl_input_required", message="Mol2MolのRLには入力SMILESが必要です。")
    return values, results, compatibility_rejections


def create_rl_job(settings: Settings, request: RlJobCreateRequest) -> JobRecord:
    if not available_torch_devices(settings).get(request.device, False):
        raise AppError(status_code=422, code="rl_device_unavailable", message=f"{request.device.upper()}は現在のREINVENT4実行環境で利用できません。")
    prior_key, model_type, prior_path, agent_path = _generation_model(settings, request.generation_model_id)
    smiles, input_results, compatibility_rejections = _input_smiles(settings, request, model_type, agent_path)
    qsar_repository = QsarModelRepository(settings)
    qsar_components = []
    for component in request.qsar_components:
        try:
            model = qsar_repository.get(component.model_id)
        except KeyError as error:
            raise AppError(status_code=404, code="qsar_model_not_found", message="選択したQSARモデルが見つかりません。") from error
        qsar_components.append({
            **component.model_dump(mode="json"),
            "model_name": model.name,
            "model_path": model.file_path,
            "task_type": model.task_type,
            "ad_threshold": model.ad_threshold,
            "algorithm": model.algorithm,
            "dataset_path": model.dataset_path,
            "evaluation_path": model.evaluation_path,
            "external_runtime_python": (
                str(settings.chemprop_runtime_python) if model.algorithm == "chemprop_v2"
                else str(settings.qsartuna_runtime_python) if model.algorithm == "qsartuna"
                else None
            ),
        })

    settings.model_store_dir.mkdir(parents=True, exist_ok=True)
    model_name = _next_rl_model_name(settings, prior_key, request.comment)
    model_id = uuid4()
    job_id = uuid4()
    job_dir = (settings.jobs_dir / str(job_id)).resolve()
    for name in ("input", "config", "logs", "output"):
        (job_dir / name).mkdir(parents=True, exist_ok=False)
    input_path = job_dir / "input" / "rl-input.smi"
    input_path.write_text("\n".join(smiles) + ("\n" if smiles else ""), encoding="utf-8")
    accepted_path = job_dir / "input" / "standardization-accepted.csv"
    rejected_path = job_dir / "input" / "standardization-rejected.csv"
    compatible = set(smiles)
    incompatibilities = dict(compatibility_rejections)
    audit_seen: set[str] = set()
    accepted_rows = []
    rejected_rows = []
    for row in input_results:
        if not row.accepted or not row.standardized_smiles:
            rejected_rows.append((row, row.reason or "invalid_smiles"))
        elif row.standardized_smiles in audit_seen:
            rejected_rows.append((row, "duplicate"))
        elif row.standardized_smiles not in compatible:
            rejected_rows.append((row, f"model_incompatible: {incompatibilities.get(row.standardized_smiles, 'unknown')}"))
            audit_seen.add(row.standardized_smiles)
        else:
            accepted_rows.append(row)
            audit_seen.add(row.standardized_smiles)
    with accepted_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Original_SMILES", "Standardized_SMILES", "MolecularWeight", "FormalCharge", "History"])
        writer.writerows((row.original_smiles, row.standardized_smiles, row.molecular_weight, row.formal_charge, ";".join(row.history)) for row in accepted_rows)
    with rejected_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Original_SMILES", "Standardized_SMILES", "MolecularWeight", "FormalCharge", "Reason", "History"])
        writer.writerows((row.original_smiles, row.standardized_smiles or "", row.molecular_weight, row.formal_charge, reason, ";".join(row.history)) for row, reason in rejected_rows)
    output_model = settings.model_store_dir / model_name
    output_model.touch(exist_ok=False)
    conditions = {
        **request.model_dump(mode="json"),
        "source": "reinforcement_learning",
        "model_name": model_name,
        "model_type": model_type,
        "base_model_key": prior_key,
        "prior_path": prior_path,
        "agent_path": agent_path,
        "agent_model_id": request.generation_model_id,
        "qsar_components": qsar_components,
        "input_smiles_path": relative_to_repository(settings, input_path),
        "output_model_path": relative_to_repository(settings, output_model),
        "summary_csv_path": relative_to_repository(settings, job_dir / "output" / "rl-summary_1.csv"),
        "smiles_standardization": standardization_metadata(),
        "standardization_counts": {
            **standardization_counts(input_results),
            "duplicate_count": sum(row.accepted for row in input_results) - len({row.standardized_smiles for row in input_results if row.accepted}),
            "model_compatibility_error_count": len(compatibility_rejections),
            "final_accepted_count": len(smiles),
        },
        "standardization_accepted_path": relative_to_repository(settings, accepted_path),
        "standardization_rejected_path": relative_to_repository(settings, rejected_path),
    }
    snapshot = job_dir / "config" / "rl-request.json"
    snapshot.write_text(json.dumps(conditions, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    conditions["request_snapshot_path"] = relative_to_repository(settings, snapshot)
    return JobRepository(settings).create(
        job_id=job_id,
        comment_original=request.comment,
        comment_slug=request.comment,
        parent_smiles="",
        run_conditions=conditions,
        job_dir=relative_to_repository(settings, job_dir),
        job_type="reinforcement_learning",
        model_id=model_id,
    )


def predict_qsar(settings: Settings, request: QsarPredictionRequest) -> list[dict]:
    try:
        model = QsarModelRepository(settings).get(request.model_id)
    except KeyError as error:
        raise AppError(status_code=404, code="qsar_model_not_found", message="QSARモデルが見つかりません。") from error
    model_path = (settings.repository_root / model.file_path).resolve()
    if not model_path.is_relative_to(settings.qsar_model_store_dir) or not model_path.is_file():
        raise AppError(status_code=500, code="invalid_qsar_model_path", message="QSARモデルファイルが不正です。")
    external = model.algorithm in {"chemprop_v2", "qsartuna"}
    script = Path(__file__).resolve().parents[2] / "worker" / ("external_qsar_predictor.py" if external else "qsar_predictor.py")
    runtime = (
        settings.chemprop_runtime_python if model.algorithm == "chemprop_v2"
        else settings.qsartuna_runtime_python if model.algorithm == "qsartuna"
        else settings.reinvent_runtime_python
    )
    argv = [str(runtime), str(script), "--model", str(model_path), "--ad-threshold", str(model.ad_threshold)]
    input_text = json.dumps({"smiles": request.smiles})
    if external:
        dataset_path = (settings.repository_root / model.dataset_path).resolve()
        metadata_path = (settings.repository_root / model.evaluation_path).resolve()
        argv.extend(["--algorithm", model.algorithm, "--training-dataset", str(dataset_path), "--metadata", str(metadata_path)])
        input_text = "\n".join(request.smiles) + "\n"
    environment = {
        "PYTHONPATH": str(Path(__file__).resolve().parents[2]),
        "PATH": os.environ.get("PATH", ""),
        "LANG": os.environ.get("LANG", "C.UTF-8"),
    }
    try:
        completed = subprocess.run(
            argv,
            input=input_text,
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
            env=environment,
        )
        payload = json.loads(completed.stdout)
        values = payload["payload"]["predictions"] if external else payload["predictions"]
        if external:
            raw_values = payload["payload"].get("raw_predictions", values)
            in_domain = payload["payload"].get("in_ad", [True] * len(values))
            return [{
                "smiles": value,
                "standardized_smiles": payload["payload"].get("standardized_smiles", [None] * len(values))[index],
                "molecular_weight": payload["payload"].get("molecular_weight", [None] * len(values))[index],
                "formal_charge": payload["payload"].get("formal_charge", [None] * len(values))[index],
                "valid": raw is not None,
                "prediction": raw,
                "scoring_prediction": score,
                "nearest_tanimoto": None,
                "in_ad": allowed,
                "error": payload["payload"].get("exclusion_reason", [None] * len(values))[index] or (None if raw is not None else "non_finite_prediction"),
            } for index, (value, raw, score, allowed) in enumerate(zip(request.smiles, raw_values, values, in_domain))]
        return values
    except (OSError, subprocess.SubprocessError, KeyError, json.JSONDecodeError) as error:
        raise AppError(status_code=500, code="qsar_prediction_failed", message="QSAR予測を実行できませんでした。") from error


def _score_transform(value: float, transform: str, low: float, high: float, k: float) -> float:
    if not math.isfinite(value):
        return 0.0
    if transform == "none":
        return value
    if transform == "double_sigmoid":
        center = (low + high) / 2
        scale = 20 / max(abs(high - low), 1.0)
        if value < center:
            exponent = max(-300.0, min(300.0, -scale * (value - low)))
            return 1 / (1 + 10 ** exponent)
        exponent = max(-300.0, min(300.0, -scale * (value - high)))
        return 1 - 1 / (1 + 10 ** exponent)
    scale = 10 * k / (high - low) if high != low else 10 * k
    exponent = max(-300.0, min(300.0, -scale * (value - (high + low) / 2)))
    sigmoid = 1 / (1 + 10 ** exponent)
    return 1 - sigmoid if transform == "reverse_sigmoid" else sigmoid


def _weighted_score(values: list[tuple[float, float]], aggregation: str) -> float:
    total_weight = sum(weight for _, weight in values)
    if not values or total_weight <= 0:
        return 0.0
    if aggregation == "arithmetic_mean":
        return sum(score * weight / total_weight for score, weight in values)
    return math.prod(max(score, 1e-8) ** (weight / total_weight) for score, weight in values)


def _property_value(molecule: Chem.Mol, component: str) -> float:
    functions = {
        "QED": lambda: QED.qed(molecule),
        "MolecularWeight": lambda: Descriptors.MolWt(molecule),
        "SlogP": lambda: Crippen.MolLogP(molecule),
        "TPSA": lambda: rdMolDescriptors.CalcTPSA(molecule),
        "HBondDonors": lambda: Lipinski.NumHDonors(molecule),
        "HBondAcceptors": lambda: Lipinski.NumHAcceptors(molecule),
        "NumRotBond": lambda: Lipinski.NumRotatableBonds(molecule),
        "SAScore": lambda: sascorer.calculateScore(molecule),
    }
    if component not in functions:
        raise AppError(status_code=422, code="rl_preview_component_unsupported", message=f"{component}は事前検証に対応していません。")
    return float(functions[component]())


def preview_rl_scoring(settings: Settings, request: RlPreviewRequest) -> list[dict]:
    qsar_repository = QsarModelRepository(settings)
    qsar_values: list[tuple[object, str, list[dict]]] = []
    for component in request.qsar_components:
        try:
            model = qsar_repository.get(component.model_id)
        except KeyError as error:
            raise AppError(status_code=404, code="qsar_model_not_found", message="QSARモデルが見つかりません。") from error
        predictions = predict_qsar(settings, QsarPredictionRequest(model_id=component.model_id, smiles=request.smiles))
        qsar_values.append((component, model.target_name or model.name, predictions))

    rows: list[dict] = []
    for index, smiles in enumerate(request.smiles):
        standardized = standardize_smiles(smiles)
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(standardized.standardized_smiles) if standardized.accepted and standardized.standardized_smiles else None
        components: list[dict] = []
        weighted: list[tuple[float, float]] = []
        for component, name, predictions in qsar_values:
            prediction = predictions[index].get("prediction")
            in_ad = predictions[index].get("in_ad") is not False
            raw = float(prediction) if isinstance(prediction, (int, float)) else math.nan
            score = 0.0 if not standardized.accepted or (not in_ad and component.ad_behavior == "penalize") else _score_transform(raw, component.transform, component.low, component.high, component.k)
            components.append({"kind": "qsar", "name": name, "prediction": prediction if math.isfinite(raw) else None, "in_ad": in_ad, "score": score, "weight": component.weight})
            weighted.append((score, component.weight))
        for component in request.property_components:
            raw = _property_value(molecule, component.component) if molecule is not None else math.nan
            score = 0.0 if not standardized.accepted else _score_transform(raw, component.transform, component.low, component.high, component.k)
            components.append({"kind": "property", "name": component.component, "prediction": raw if math.isfinite(raw) else None, "in_ad": None, "score": score, "weight": component.weight})
            weighted.append((score, component.weight))
        rows.append({
            "smiles": smiles,
            "generated_smiles": smiles,
            "standardized_smiles": standardized.standardized_smiles,
            "molecular_weight": standardized.molecular_weight,
            "formal_charge": standardized.formal_charge,
            "valid": standardized.accepted,
            "exclusion_reason": standardized.reason,
            "components": components,
            "total_score": _weighted_score(weighted, request.aggregation) if standardized.accepted else 0.0,
        })
    return rows
