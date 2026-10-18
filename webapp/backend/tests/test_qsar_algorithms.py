from __future__ import annotations

import importlib.util

import numpy as np
import pytest

pytest.importorskip("joblib")
pytest.importorskip("sklearn")

from worker.qsar_training import _pipeline


REGRESSION_SMILES = np.asarray([
    "CC", "CCC", "CCCC", "CCO", "CCCO", "CCN", "CCCN", "c1ccccc1",
    "Cc1ccccc1", "Oc1ccccc1", "C1CCCCC1", "C1CCCC1",
], dtype=object)
REGRESSION_TARGETS = np.linspace(0.1, 1.2, len(REGRESSION_SMILES))
CLASSIFICATION_TARGETS = np.asarray([0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1])


def config(algorithm: str, task_type: str) -> dict:
    return {
        "algorithm": algorithm,
        "task_type": task_type,
        "fingerprint": "morgan",
        "feature_mode": "fingerprint",
        "morgan_radius": 2,
        "bit_length": 128,
        "use_chirality": True,
        "advanced": False,
        "hyperparameters": {},
    }


@pytest.mark.parametrize("algorithm", ["kernel_ridge", "gaussian_process", "ridge", "elastic_net", "mlp"])
def test_added_sklearn_regressors_fit_predict(algorithm):
    pipeline = _pipeline(config(algorithm, "regression"))
    pipeline.fit(REGRESSION_SMILES, REGRESSION_TARGETS)
    prediction = pipeline.predict(REGRESSION_SMILES[:2])
    assert prediction.shape == (2,)
    assert np.isfinite(prediction).all()


@pytest.mark.parametrize("algorithm", ["gaussian_process", "logistic_regression", "mlp"])
def test_added_sklearn_classifiers_fit_predict_proba(algorithm):
    pipeline = _pipeline(config(algorithm, "classification"))
    pipeline.fit(REGRESSION_SMILES, CLASSIFICATION_TARGETS)
    probabilities = pipeline.predict_proba(REGRESSION_SMILES[:2])
    assert probabilities.shape == (2, 2)
    assert np.isfinite(probabilities).all()


@pytest.mark.parametrize("algorithm,package", [
    ("lightgbm", "lightgbm"),
    ("xgboost", "xgboost"),
    ("catboost", "catboost"),
])
@pytest.mark.parametrize("task_type", ["regression", "classification"])
def test_external_boosters_fit_and_reload_compatible(algorithm, package, task_type):
    if importlib.util.find_spec(package) is None:
        pytest.skip(f"{package} is not installed")
    targets = REGRESSION_TARGETS if task_type == "regression" else CLASSIFICATION_TARGETS
    pipeline = _pipeline(config(algorithm, task_type))
    pipeline.set_params(model__n_estimators=5)
    pipeline.fit(REGRESSION_SMILES, targets)
    prediction = pipeline.predict(REGRESSION_SMILES[:2])
    assert prediction.shape == (2,)
    if task_type == "classification":
        assert pipeline.predict_proba(REGRESSION_SMILES[:2]).shape == (2, 2)
