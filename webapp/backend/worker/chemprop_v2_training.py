from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
from rdkit import DataStructs
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, mean_absolute_error, mean_squared_error, precision_score, r2_score, recall_score, roc_auc_score
from sklearn.model_selection import train_test_split

from worker.qsar_features import ad_fingerprint
from worker.qsar_training import _read_dataset, _scaffold_split


def _indices(smiles: np.ndarray, targets: np.ndarray, config: dict):
    if config["split_method"] == "scaffold":
        train_val, test = _scaffold_split(smiles, config["train_fraction"])
        inner_train, val = _scaffold_split(smiles[train_val], 0.8)
        return train_val[inner_train], train_val[val], test
    indices = np.arange(len(smiles))
    stratify = targets if config["task_type"] == "classification" else None
    train_val, test = train_test_split(indices, train_size=config["train_fraction"], random_state=42, stratify=stratify)
    inner_stratify = targets[train_val] if stratify is not None else None
    train, val = train_test_split(train_val, train_size=0.8, random_state=42, stratify=inner_stratify)
    return train, val, test


def _metrics(actual, predicted, task_type: str):
    if task_type == "regression":
        return {"r2": float(r2_score(actual, predicted)) if len(actual) > 1 else None, "rmse": float(math.sqrt(mean_squared_error(actual, predicted))), "mae": float(mean_absolute_error(actual, predicted))}
    labels = (np.asarray(predicted, dtype=float) >= 0.5).astype(int)
    actual = np.asarray(actual, dtype=int)
    result = {
        "accuracy": float(accuracy_score(actual, labels)),
        "precision": float(precision_score(actual, labels, zero_division=0)),
        "recall": float(recall_score(actual, labels, zero_division=0)),
        "f1": float(f1_score(actual, labels, zero_division=0)),
    }
    if len(np.unique(actual)) == 2:
        result.update({"roc_auc": float(roc_auc_score(actual, predicted)), "pr_auc": float(average_precision_score(actual, predicted))})
    else:
        result.update({"roc_auc": None, "pr_auc": None})
    return result


def train(config_path: Path, repository_root: Path) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    runtime = Path(config["external_runtime_python"])
    chemprop = runtime.with_name("chemprop")
    if not chemprop.is_file():
        raise RuntimeError(f"Chemprop executable was not found next to {runtime}")
    dataset_path = repository_root / config["dataset_path"]
    output_model = repository_root / config["output_model_path"]
    predictions_path = repository_root / config["predictions_path"]
    evaluation_path = repository_root / config["evaluation_path"]
    work_dir = config_path.parent.parent / "output" / "chemprop-v2"
    work_dir.mkdir(parents=True, exist_ok=True)
    smiles, targets = _read_dataset(dataset_path, config["task_type"])
    class_mapping = None
    if config["task_type"] == "classification":
        classes = sorted(set(str(value) for value in targets))
        if len(classes) != 2:
            raise ValueError("Chemprop classification currently requires exactly two classes")
        class_mapping = {value: index for index, value in enumerate(classes)}
        targets = np.asarray([class_mapping[str(value)] for value in targets], dtype=int)
    train_idx, val_idx, test_idx = _indices(smiles, targets, config)
    prepared = work_dir / "dataset.csv"
    with prepared.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["smiles", "Target", "split"])
        labels = np.full(len(smiles), "train", dtype=object)
        labels[val_idx] = "val"
        labels[test_idx] = "test"
        writer.writerows(zip(smiles, targets, labels))
    hp = config.get("hyperparameters", {}) if config.get("advanced") else {}
    epochs = int(hp.get("epochs", 30))
    accelerator = {"cuda": "gpu", "mps": "mps"}.get(config.get("device"), "cpu")
    command = [
        str(chemprop), "train", "--data-path", str(prepared),
        "--smiles-columns", "smiles", "--target-columns", "Target", "--splits-column", "split",
        "--task-type", config["task_type"], "--output-dir", str(work_dir / "checkpoints"),
        "--epochs", str(epochs), "--warmup-epochs", str(min(2, max(0, epochs - 1))),
        "--batch-size", str(int(hp.get("batch_size", 64))),
        "--init-lr", str(float(hp.get("learning_rate", 1e-4))), "--patience", str(int(hp.get("patience", 10))),
        "--accelerator", accelerator, "--devices", "1", "--num-workers", "0", "--save-data-splits",
    ]
    subprocess.run(command, check=True)
    checkpoints = sorted((work_dir / "checkpoints").rglob("*.ckpt"))
    if not checkpoints:
        raise RuntimeError("Chemprop completed without a checkpoint")
    checkpoint = next((path for path in checkpoints if "best" in path.name.lower()), checkpoints[0])
    output_model.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(checkpoint, output_model)

    test_csv = work_dir / "test.csv"
    with test_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["smiles"])
        writer.writerows((smiles[index],) for index in test_idx)
    raw_predictions = work_dir / "raw-predictions.csv"
    subprocess.run([str(chemprop), "predict", "--test-path", str(test_csv), "--smiles-columns", "smiles", "--model-path", str(output_model), "--preds-path", str(raw_predictions), "--accelerator", accelerator, "--devices", "1", "--num-workers", "0"], check=True)
    with raw_predictions.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    numeric_columns = [column for column in (rows[0].keys() if rows else []) if column != "smiles"]
    if not numeric_columns:
        raise RuntimeError("Chemprop prediction output has no prediction column")
    predicted = np.asarray([float(row[numeric_columns[0]]) for row in rows], dtype=float)
    actual = targets[test_idx]
    train_fps = [ad_fingerprint(str(smiles[index])) for index in train_idx]
    similarities = []
    for index in test_idx:
        fingerprint = ad_fingerprint(str(smiles[index]))
        scores = DataStructs.BulkTanimotoSimilarity(fingerprint, train_fps) if fingerprint else []
        similarities.append(max(scores) if scores else 0.0)
    predictions_path.parent.mkdir(parents=True, exist_ok=True)
    with predictions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Split", "SMILES", "Actual", "Predicted", "NearestTanimoto", "In_AD"])
        for index, prediction, similarity in zip(test_idx, predicted, similarities):
            writer.writerow(["test", smiles[index], targets[index], prediction, similarity, similarity >= config["ad_threshold"]])
    test_metrics = _metrics(actual, predicted, config["task_type"])
    metrics = {
        "training": {}, "test": test_metrics,
        "cross_validation": {"available": False, "reason": "Chemprop uses an explicit validation split and early stopping."},
        "training_data_count": len(train_idx), "validation_data_count": len(val_idx), "test_data_count": len(test_idx),
    }
    plots = {"actual_vs_predicted": [[float(a), float(p)] for a, p in zip(actual, predicted)]} if config["task_type"] == "regression" else {}
    metadata = {
        **config, "algorithm": "chemprop_v2", "chemprop_version": subprocess.run([str(chemprop), "-v"], check=True, capture_output=True, text=True).stdout.strip().split()[-1],
        "class_mapping": class_mapping, "training_smiles": [str(smiles[index]) for index in train_idx],
        "ad_method": "Morgan radius 2 / 2048-bit nearest-neighbour Tanimoto",
    }
    evaluation_path.write_text(json.dumps({"metrics": metrics, "metadata": metadata, "plots": plots}, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    args = parser.parse_args()
    train(args.config, args.repository_root.resolve())


if __name__ == "__main__":
    main()
