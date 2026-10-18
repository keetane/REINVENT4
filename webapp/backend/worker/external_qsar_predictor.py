from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from app.services.smiles_standardization import assert_standardization_compatible, standardize_smiles


def _ad_scores(smiles: list[str], dataset: Path, threshold: float) -> list[bool]:
    if threshold <= 0:
        return [True] * len(smiles)
    with dataset.open(encoding="utf-8", newline="") as handle:
        training = [row["Canonical_SMILES"] for row in csv.DictReader(handle)]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    train_fps = [generator.GetFingerprint(molecule) for value in training if (molecule := Chem.MolFromSmiles(value)) is not None]
    result = []
    for value in smiles:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            result.append(False)
            continue
        scores = DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(molecule), train_fps)
        result.append(bool(scores and max(scores) >= threshold))
    return result


def _chemprop(model: Path, smiles: list[str]) -> list[float]:
    with tempfile.TemporaryDirectory(prefix="chemprop-predict-") as directory:
        root = Path(directory)
        data = root / "input.csv"
        output = root / "predictions.csv"
        with data.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["smiles"])
            writer.writerows((value,) for value in smiles)
        chemprop = Path(sys.executable).with_name("chemprop")
        if not chemprop.is_file():
            raise RuntimeError(f"Chemprop executable was not found next to {sys.executable}")
        subprocess.run([str(chemprop), "predict", "--test-path", str(data), "--smiles-columns", "smiles", "--model-path", str(model), "--preds-path", str(output), "--accelerator", "cpu", "--devices", "1", "--num-workers", "0", "-qq"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with output.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        columns = [column for column in rows[0] if column != "smiles"] if rows else []
        if not columns:
            raise RuntimeError("Chemprop produced no prediction column")
        return [float(row[columns[0]]) for row in rows]


def _qsartuna(model: Path, smiles: list[str]) -> list[float]:
    import dill

    with model.open("rb") as handle:
        predictor = dill.load(handle)
    return [float(value) for value in predictor.predict_from_smiles(smiles, uncert=False)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--algorithm", choices=["chemprop_v2", "qsartuna"], required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--training-dataset", type=Path, required=True)
    parser.add_argument("--ad-threshold", type=float, default=0.0)
    parser.add_argument("--below-ad-score", type=float, default=0.0)
    parser.add_argument("--metadata", type=Path, required=True)
    args = parser.parse_args()
    smiles = [line.strip() for line in sys.stdin if line.strip()]
    evaluation = json.loads(args.metadata.read_text(encoding="utf-8"))
    assert_standardization_compatible(evaluation.get("metadata"))
    standardized = [standardize_smiles(value) for value in smiles]
    valid_indices = [index for index, result in enumerate(standardized) if result.accepted and result.standardized_smiles]
    valid_smiles = [standardized[index].standardized_smiles for index in valid_indices]
    predictions = (
        _chemprop(args.model, valid_smiles) if args.algorithm == "chemprop_v2" else _qsartuna(args.model, valid_smiles)
    ) if valid_smiles else []
    valid_in_domain = _ad_scores(valid_smiles, args.training_dataset, args.ad_threshold)
    finite = [math.isfinite(float(value)) for value in predictions]
    raw_predictions: list[float | None] = [None] * len(smiles)
    safe: list[float] = [0.0] * len(smiles)
    in_domain: list[bool] = [False] * len(smiles)
    for source_index, value, allowed, valid in zip(valid_indices, predictions, valid_in_domain, finite):
        raw_predictions[source_index] = float(value) if valid else None
        safe[source_index] = float(value) if allowed and valid else args.below_ad_score
        in_domain[source_index] = allowed and valid
    print(json.dumps({"version": 1, "payload": {
        "predictions": safe,
        "raw_predictions": raw_predictions,
        "in_ad": in_domain,
        "standardized_smiles": [result.standardized_smiles for result in standardized],
        "molecular_weight": [result.molecular_weight for result in standardized],
        "formal_charge": [result.formal_charge for result in standardized],
        "exclusion_reason": [result.reason for result in standardized],
    }}))


if __name__ == "__main__":
    main()
