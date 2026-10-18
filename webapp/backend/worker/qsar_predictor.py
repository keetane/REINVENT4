from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import joblib
import numpy as np
from rdkit import DataStructs

from app.services.smiles_standardization import assert_standardization_compatible, standardize_smiles
from worker.qsar_features import ad_fingerprint


def predict(model_path: Path, smiles: list[str], ad_threshold: float):
    bundle = joblib.load(model_path)
    pipeline = bundle["pipeline"]
    metadata = bundle["metadata"]
    assert_standardization_compatible(metadata)
    training_fps = [ad_fingerprint(value) for value in bundle["training_smiles"]]
    valid_indices, valid_smiles = [], []
    standardized = [standardize_smiles(value) for value in smiles]
    rows = [{
        "smiles": value,
        "standardized_smiles": result.standardized_smiles,
        "molecular_weight": result.molecular_weight,
        "formal_charge": result.formal_charge,
        "valid": False,
        "prediction": None,
        "scoring_prediction": 0.0,
        "nearest_tanimoto": None,
        "in_ad": False,
        "error": result.reason,
    } for value, result in zip(smiles, standardized)]
    for index, value in enumerate(smiles):
        result = standardized[index]
        if result.accepted and result.standardized_smiles:
            valid_indices.append(index)
            valid_smiles.append(result.standardized_smiles)
    if not valid_smiles:
        return rows
    values = np.asarray(valid_smiles, dtype=object)
    if metadata["task_type"] == "classification":
        probabilities = pipeline.predict_proba(values)
        predictions = probabilities[:, 1] if probabilities.shape[1] == 2 else probabilities.max(axis=1)
    else:
        predictions = pipeline.predict(values)
    for source_index, value, prediction in zip(valid_indices, valid_smiles, predictions):
        fingerprint = ad_fingerprint(value)
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprint, training_fps) if fingerprint else []
        nearest = max(similarities) if similarities else 0.0
        finite = bool(np.isfinite(prediction))
        rows[source_index] = {
            "smiles": smiles[source_index],
            "standardized_smiles": value,
            "molecular_weight": standardized[source_index].molecular_weight,
            "formal_charge": standardized[source_index].formal_charge,
            "valid": finite,
            "prediction": float(prediction) if finite else None,
            "scoring_prediction": float(prediction) if finite else 0.0,
            "nearest_tanimoto": float(nearest),
            "in_ad": finite and nearest >= ad_threshold,
            "error": None if finite else "non_finite_prediction",
        }
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--ad-threshold", required=True, type=float)
    arguments = parser.parse_args()
    payload = json.load(sys.stdin)
    print(json.dumps({"predictions": predict(arguments.model, payload["smiles"], arguments.ad_threshold)}))


if __name__ == "__main__":
    main()
