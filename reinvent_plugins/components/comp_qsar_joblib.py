"""Trusted scikit-learn QSAR models as a REINVENT4 scoring component."""

from __future__ import annotations

__all__ = ["QsarJoblib"]

from typing import List

import joblib
import numpy as np
from pydantic.dataclasses import dataclass
from rdkit import DataStructs

from worker.qsar_features import ad_fingerprint
from app.services.smiles_standardization import assert_standardization_compatible, standardize_smiles

from .add_tag import add_tag
from .component_results import ComponentResults


@add_tag("__parameters")
@dataclass
class Parameters:
    model_file: List[str]
    task_type: List[str]
    ad_threshold: List[float]
    below_ad_score: List[float]
    class_index: List[int]


@add_tag("__component")
class QsarJoblib:
    """Batch predictions from server-created, trusted joblib bundles."""

    def __init__(self, params: Parameters):
        self.models = []
        for index, filename in enumerate(params.model_file):
            bundle = joblib.load(filename)
            if not isinstance(bundle, dict) or "pipeline" not in bundle or "training_smiles" not in bundle:
                raise ValueError("invalid QSAR model bundle")
            assert_standardization_compatible(bundle.get("metadata"))
            training_fps = [ad_fingerprint(value) for value in bundle["training_smiles"]]
            self.models.append(
                (
                    bundle["pipeline"],
                    training_fps,
                    params.task_type[index],
                    float(params.ad_threshold[index]),
                    float(params.below_ad_score[index]),
                    int(params.class_index[index]),
                )
            )
        self.number_of_endpoints = len(self.models)

    def __call__(self, smilies: List[str]) -> ComponentResults:
        endpoint_scores = []
        for pipeline, training_fps, task_type, threshold, below_ad_score, class_index in self.models:
            standardized = [standardize_smiles(value) for value in smilies]
            valid_indices = [index for index, result in enumerate(standardized) if result.accepted and result.standardized_smiles]
            valid_smiles = [standardized[index].standardized_smiles for index in valid_indices]
            scores = np.zeros(len(smilies), dtype=float)
            if valid_smiles and task_type == "classification":
                probabilities = pipeline.predict_proba(np.asarray(valid_smiles, dtype=object))
                if class_index < 0 or class_index >= probabilities.shape[1]:
                    raise ValueError("classification class index is outside model output")
                predicted = np.asarray(probabilities[:, class_index], dtype=float)
            elif valid_smiles:
                predicted = np.asarray(pipeline.predict(np.asarray(valid_smiles, dtype=object)), dtype=float)
            else:
                predicted = np.asarray([], dtype=float)
            for row, smiles, prediction in zip(valid_indices, valid_smiles, predicted):
                scores[row] = prediction
                fingerprint = ad_fingerprint(smiles)
                if fingerprint is None:
                    scores[row] = 0.0
                    continue
                similarities = DataStructs.BulkTanimotoSimilarity(fingerprint, training_fps)
                if not similarities or max(similarities) < threshold:
                    scores[row] = below_ad_score
            scores[~np.isfinite(scores)] = 0.0
            endpoint_scores.append(scores)
        return ComponentResults(endpoint_scores)
