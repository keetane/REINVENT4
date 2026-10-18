"""Applicability-domain penalty evaluated on standardized structures."""

from __future__ import annotations

__all__ = ["QsarApplicabilityDomain"]

import csv
from pathlib import Path
from typing import List, Optional

import joblib
import numpy as np
from pydantic.dataclasses import dataclass
from rdkit import DataStructs

from app.services.smiles_standardization import assert_standardization_compatible, standardize_smiles
from worker.qsar_features import ad_fingerprint

from .add_tag import add_tag
from .component_results import ComponentResults


@add_tag("__parameters")
@dataclass
class Parameters:
    model_file: List[Optional[str]]
    dataset_file: List[Optional[str]]
    threshold: List[float]


@add_tag("__component", "penalty")
class QsarApplicabilityDomain:
    def __init__(self, params: Parameters):
        self.endpoints = []
        for model_file, dataset_file, threshold in zip(params.model_file, params.dataset_file, params.threshold):
            if model_file:
                bundle = joblib.load(model_file)
                assert_standardization_compatible(bundle.get("metadata"))
                training_smiles = bundle["training_smiles"]
            elif dataset_file:
                with Path(dataset_file).open(encoding="utf-8", newline="") as handle:
                    training_smiles = [row["Canonical_SMILES"] for row in csv.DictReader(handle)]
            else:
                raise ValueError("QSAR AD requires a model or training dataset")
            self.endpoints.append(([ad_fingerprint(value) for value in training_smiles], float(threshold)))
        self.number_of_endpoints = len(self.endpoints)

    def __call__(self, smilies: List[str]) -> ComponentResults:
        endpoint_scores = []
        for training_fps, threshold in self.endpoints:
            scores = []
            for value in smilies:
                result = standardize_smiles(value)
                fingerprint = ad_fingerprint(result.standardized_smiles) if result.accepted and result.standardized_smiles else None
                similarities = DataStructs.BulkTanimotoSimilarity(fingerprint, training_fps) if fingerprint else []
                scores.append(1.0 if similarities and max(similarities) >= threshold else 0.0)
            endpoint_scores.append(np.asarray(scores, dtype=float))
        return ComponentResults(endpoint_scores)
