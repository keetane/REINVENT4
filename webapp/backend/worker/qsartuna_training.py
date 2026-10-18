from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
from rdkit import DataStructs
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, mean_absolute_error, mean_squared_error, precision_score, r2_score, recall_score, roc_auc_score

from worker.qsar_features import ad_fingerprint
from worker.qsar_training import _read_dataset, _split


def _metrics(actual, predicted, task_type):
    if task_type == "regression":
        return {"r2": float(r2_score(actual, predicted)) if len(actual) > 1 else None, "rmse": float(math.sqrt(mean_squared_error(actual, predicted))), "mae": float(mean_absolute_error(actual, predicted))}
    labels = (np.asarray(predicted) >= 0.5).astype(int)
    actual = np.asarray(actual, dtype=int)
    result = {"accuracy": float(accuracy_score(actual, labels)), "precision": float(precision_score(actual, labels, zero_division=0)), "recall": float(recall_score(actual, labels, zero_division=0)), "f1": float(f1_score(actual, labels, zero_division=0))}
    if len(np.unique(actual)) == 2:
        result.update({"roc_auc": float(roc_auc_score(actual, predicted)), "pr_auc": float(average_precision_score(actual, predicted))})
    return result


def train(config_path: Path, repository_root: Path) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    runtime = Path(config["external_runtime_python"])
    dataset_path = repository_root / config["dataset_path"]
    output_model = repository_root / config["output_model_path"]
    predictions_path = repository_root / config["predictions_path"]
    evaluation_path = repository_root / config["evaluation_path"]
    smiles, targets = _read_dataset(dataset_path, config["task_type"])
    class_mapping = None
    if config["task_type"] == "classification":
        classes = sorted(set(str(value) for value in targets))
        if len(classes) != 2:
            raise ValueError("QSARtuna classification currently requires exactly two classes")
        class_mapping = {value: index for index, value in enumerate(classes)}
        targets = np.asarray([class_mapping[str(value)] for value in targets], dtype=int)
    train_smiles, test_smiles, train_targets, test_targets = _split(smiles, targets, config)
    output_model.parent.mkdir(parents=True, exist_ok=True)
    work_dir = config_path.parent.parent / "output" / "qsartuna"
    work_dir.mkdir(parents=True, exist_ok=True)
    train_csv, test_csv = work_dir / "train.csv", work_dir / "test.csv"
    for path, values, labels in ((train_csv, train_smiles, train_targets), (test_csv, test_smiles, test_targets)):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle); writer.writerow(["Canonical_SMILES", "Target"]); writer.writerows(zip(values, labels))
    hp = config.get("hyperparameters", {}) if config.get("advanced") else {}
    timeout_seconds = int(hp.get("timeout_seconds", 3600))
    runner = Path(__file__).with_name("qsartuna_runner.py")
    build_config = work_dir / "best-build-config.json"
    environment = {**os.environ, "NUMBA_DISABLE_JIT": "1", "MPLCONFIGDIR": str(work_dir / "matplotlib")}
    subprocess.run([
        str(runtime), str(runner), "--train", str(train_csv), "--test", str(test_csv),
        "--model", str(output_model), "--build-config", str(build_config), "--task-type", config["task_type"],
        "--split-method", config["split_method"], "--folds", str(config["cross_validation_folds"]),
        "--trials", str(int(hp.get("n_trials", 30))), "--startup-trials", str(int(hp.get("n_startup_trials", 10))),
    ], check=True, env=environment, timeout=timeout_seconds)
    predictor = Path(__file__).with_name("external_qsar_predictor.py")
    predictor_metadata = work_dir / "predictor-metadata.json"
    predictor_metadata.write_text(json.dumps({"metadata": config}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run([
        str(runtime), str(predictor), "--algorithm", "qsartuna", "--model", str(output_model),
        "--training-dataset", str(train_csv), "--ad-threshold", "0", "--metadata", str(predictor_metadata),
    ], input="\n".join(str(value) for value in test_smiles) + "\n", text=True, capture_output=True, check=True, env=environment)
    predicted = np.asarray(json.loads(completed.stdout)["payload"]["predictions"], dtype=float)
    train_fps = [ad_fingerprint(str(value)) for value in train_smiles]
    similarities = []
    for value in test_smiles:
        fingerprint = ad_fingerprint(str(value)); scores = DataStructs.BulkTanimotoSimilarity(fingerprint, train_fps) if fingerprint else []
        similarities.append(max(scores) if scores else 0.0)
    with predictions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle); writer.writerow(["Split", "SMILES", "Actual", "Predicted", "NearestTanimoto", "In_AD"])
        for value, actual, prediction, similarity in zip(test_smiles, test_targets, predicted, similarities):
            writer.writerow(["test", value, actual, prediction, similarity, similarity >= config["ad_threshold"]])
    metrics = {"training": {}, "test": _metrics(test_targets, predicted, config["task_type"]), "cross_validation": {"available": True, "source": "QSARtuna Optuna study"}, "training_data_count": len(train_smiles), "test_data_count": len(test_smiles)}
    plots = {"actual_vs_predicted": [[float(a), float(p)] for a, p in zip(test_targets, predicted)]} if config["task_type"] == "regression" else {}
    version = subprocess.run([str(runtime), "-c", "import importlib.metadata as m; print(m.version('qsartuna'))"], check=True, capture_output=True, text=True, env=environment).stdout.strip()
    metadata = {**config, "algorithm": "qsartuna", "qsartuna_version": version, "class_mapping": class_mapping, "best_build_config_path": str(build_config.relative_to(repository_root)), "training_smiles": [str(value) for value in train_smiles]}
    evaluation_path.write_text(json.dumps({"metrics": metrics, "metadata": metadata, "plots": plots}, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(); parser.add_argument("--config", required=True, type=Path); parser.add_argument("--repository-root", required=True, type=Path)
    args = parser.parse_args(); train(args.config, args.repository_root.resolve())


if __name__ == "__main__":
    main()
