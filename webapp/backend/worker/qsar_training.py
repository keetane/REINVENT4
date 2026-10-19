from __future__ import annotations

import argparse
import csv
import importlib.metadata
import json
import math
import platform
import time
from pathlib import Path

import joblib
import numpy as np
import rdkit
import sklearn
from rdkit import DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import (
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    GradientBoostingClassifier,
    GradientBoostingRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
)
from sklearn.impute import SimpleImputer
from sklearn.gaussian_process import GaussianProcessClassifier, GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import ElasticNet, LogisticRegression, Ridge
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    confusion_matrix,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    precision_recall_curve,
    r2_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold, cross_validate, train_test_split
from sklearn.neural_network import MLPClassifier, MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, SVR

from worker.qsar_features import SmilesFeaturizer, ad_fingerprint


def _package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def _read_dataset(path: Path, task_type: str):
    smiles, targets = [], []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            smiles.append(row["Canonical_SMILES"])
            targets.append(row["Target"] if task_type == "classification" else float(row["Target"]))
    return np.asarray(smiles, dtype=object), np.asarray(targets)


def _scaffold_split(smiles: np.ndarray, train_fraction: float):
    groups: dict[str, list[int]] = {}
    for index, value in enumerate(smiles):
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(smiles=str(value), includeChirality=True)
        groups.setdefault(scaffold or f"acyclic:{value}", []).append(index)
    if len(groups) < 2:
        raise ValueError("Scaffold Splitには異なるScaffoldが2種類以上必要です。Random Splitを選択してください。")
    ordered = sorted(groups.values(), key=lambda group: (-len(group), group[0]))
    target_test = max(1, round(len(smiles) * (1 - train_fraction)))
    test_indices: list[int] = []
    train_indices: list[int] = []
    for group in ordered:
        if len(test_indices) < target_test and len(test_indices) + len(group) < len(smiles):
            test_indices.extend(group)
        else:
            train_indices.extend(group)
    if not train_indices or not test_indices:
        raise ValueError("Scaffold Splitを成立させられません。Random Splitを選択してください。")
    return np.asarray(train_indices), np.asarray(test_indices)


def _split(smiles, targets, config):
    if config["split_method"] == "scaffold":
        train_idx, test_idx = _scaffold_split(smiles, config["train_fraction"])
    else:
        indices = np.arange(len(smiles))
        stratify = targets if config["task_type"] == "classification" else None
        train_idx, test_idx = train_test_split(
            indices,
            train_size=config["train_fraction"],
            random_state=42,
            stratify=stratify,
        )
    return smiles[train_idx], smiles[test_idx], targets[train_idx], targets[test_idx]


def _bounded(parameters, allowed):
    return {key: value for key, value in parameters.items() if key in allowed}


def _estimator(config):
    algorithm = config["algorithm"]
    classification = config["task_type"] == "classification"
    supplied = config.get("hyperparameters", {}) if config.get("advanced") else {}
    if algorithm == "random_forest":
        params = {"n_estimators": 300, "random_state": 42, "n_jobs": -1}
        params.update(_bounded(supplied, {"n_estimators", "max_depth", "min_samples_split", "min_samples_leaf"}))
        return RandomForestClassifier(**params) if classification else RandomForestRegressor(**params)
    if algorithm == "svm":
        params = {"C": 1.0, "kernel": "rbf", "gamma": "scale"}
        params.update(_bounded(supplied, {"C", "kernel", "gamma", "epsilon"}))
        if classification:
            params.pop("epsilon", None)
            return SVC(**params, probability=True, random_state=42)
        return SVR(**params)
    if algorithm == "gradient_boosting":
        params = {"n_estimators": 200, "learning_rate": 0.05, "max_depth": 3, "random_state": 42}
        params.update(_bounded(supplied, {"n_estimators", "learning_rate", "max_depth"}))
        return GradientBoostingClassifier(**params) if classification else GradientBoostingRegressor(**params)
    if algorithm == "extra_trees":
        params = {"n_estimators": 300, "random_state": 42, "n_jobs": -1}
        params.update(_bounded(supplied, {"n_estimators", "max_depth", "min_samples_leaf"}))
        return ExtraTreesClassifier(**params) if classification else ExtraTreesRegressor(**params)
    if algorithm == "kernel_ridge":
        if classification:
            raise ValueError("Kernel Ridgeは回帰専用です。")
        params = {"alpha": 1.0, "kernel": "rbf", "gamma": None}
        params.update(_bounded(supplied, {"alpha", "kernel", "gamma", "degree"}))
        return KernelRidge(**params)
    if algorithm == "gaussian_process":
        kernel_name = str(supplied.get("kernel", "rbf"))
        if kernel_name not in {"rbf", "matern"}:
            raise ValueError("unsupported Gaussian Process kernel")
        kernel = RBF() if kernel_name == "rbf" else Matern(nu=1.5)
        params = {"kernel": kernel, "random_state": 42}
        if classification:
            return GaussianProcessClassifier(**params)
        params["normalize_y"] = True
        return GaussianProcessRegressor(**params)
    if algorithm == "ridge":
        if classification:
            raise ValueError("Ridgeは回帰専用です。分類にはLogistic Regressionを選択してください。")
        params = {"alpha": 1.0}
        params.update(_bounded(supplied, {"alpha"}))
        return Ridge(**params)
    if algorithm == "elastic_net":
        if classification:
            raise ValueError("Elastic Netは回帰専用です。分類にはLogistic Regressionを選択してください。")
        params = {"alpha": 1.0, "l1_ratio": 0.5, "max_iter": 5000, "random_state": 42}
        params.update(_bounded(supplied, {"alpha", "l1_ratio", "max_iter"}))
        return ElasticNet(**params)
    if algorithm == "logistic_regression":
        if not classification:
            raise ValueError("Logistic Regressionは分類専用です。")
        params = {"C": 1.0, "max_iter": 5000, "random_state": 42}
        params.update(_bounded(supplied, {"C", "max_iter"}))
        return LogisticRegression(**params)
    if algorithm == "mlp":
        params = {"hidden_layer_sizes": (128, 64), "max_iter": 500, "early_stopping": True, "random_state": 42}
        params.update(_bounded(supplied, {"hidden_layer_sizes", "alpha", "learning_rate_init", "max_iter"}))
        if isinstance(params.get("hidden_layer_sizes"), list):
            params["hidden_layer_sizes"] = tuple(int(value) for value in params["hidden_layer_sizes"])
        return MLPClassifier(**params) if classification else MLPRegressor(**params)
    if algorithm == "lightgbm":
        from lightgbm import LGBMClassifier, LGBMRegressor
        params = {"n_estimators": 300, "learning_rate": 0.05, "max_depth": -1, "num_leaves": 31, "random_state": 42, "n_jobs": -1, "verbosity": -1}
        if config.get("device") == "cuda":
            params["device_type"] = "cuda"
        params.update(_bounded(supplied, {"n_estimators", "learning_rate", "max_depth", "num_leaves", "min_child_samples", "subsample", "colsample_bytree"}))
        return LGBMClassifier(**params) if classification else LGBMRegressor(**params)
    if algorithm == "xgboost":
        from xgboost import XGBClassifier, XGBRegressor
        params = {"n_estimators": 300, "learning_rate": 0.05, "max_depth": 6, "random_state": 42, "n_jobs": -1}
        if config.get("device") == "cuda":
            params.update({"tree_method": "hist", "device": "cuda"})
        params.update(_bounded(supplied, {"n_estimators", "learning_rate", "max_depth", "min_child_weight", "subsample", "colsample_bytree", "reg_alpha", "reg_lambda"}))
        if classification:
            params.update({"eval_metric": "logloss"})
            return XGBClassifier(**params)
        return XGBRegressor(**params)
    if algorithm == "catboost":
        from catboost import CatBoostClassifier, CatBoostRegressor
        params = {"iterations": 300, "learning_rate": 0.05, "depth": 6, "random_seed": 42, "verbose": False, "allow_writing_files": False}
        if config.get("device") == "cuda":
            params.update({"task_type": "GPU", "devices": "0"})
        params.update(_bounded(supplied, {"iterations", "learning_rate", "depth", "l2_leaf_reg"}))
        return CatBoostClassifier(**params) if classification else CatBoostRegressor(**params)
    raise ValueError(f"unsupported algorithm: {algorithm}")


def _pipeline(config):
    steps = [
        ("features", SmilesFeaturizer(
            fingerprint=config["fingerprint"],
            feature_mode=config["feature_mode"],
            morgan_radius=config["morgan_radius"],
            bit_length=config["bit_length"],
            use_chirality=config["use_chirality"],
        )),
        ("impute", SimpleImputer(strategy="median")),
    ]
    if config["algorithm"] in {"svm", "kernel_ridge", "gaussian_process", "ridge", "elastic_net", "logistic_regression", "mlp"}:
        steps.append(("scale", StandardScaler()))
    steps.append(("model", _estimator(config)))
    return Pipeline(steps)


def _cross_validation(pipeline, smiles, targets, config):
    folds = min(config["cross_validation_folds"], len(smiles))
    if config["task_type"] == "classification":
        classes, counts = np.unique(targets, return_counts=True)
        folds = min(folds, int(counts.min()))
        if folds < 2:
            return {"available": False, "reason": "各クラスのTrainingデータが2件未満です。"}
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=42)
        scoring = {"accuracy": "accuracy", "f1": "f1_weighted"}
        if len(classes) == 2:
            scoring.update({"roc_auc": "roc_auc", "pr_auc": "average_precision"})
    else:
        if folds < 2:
            return {"available": False, "reason": "Trainingデータが不足しています。"}
        splitter = KFold(n_splits=folds, shuffle=True, random_state=42)
        scoring = {"r2": "r2", "mae": "neg_mean_absolute_error", "rmse": "neg_root_mean_squared_error"}
    results = cross_validate(pipeline, smiles, targets, cv=splitter, scoring=scoring, error_score="raise")
    summary = {"available": True, "folds": folds}
    for key in scoring:
        values = results[f"test_{key}"]
        if key in {"mae", "rmse"}:
            values = -values
        mean = float(np.mean(values))
        summary[key] = mean
        summary[f"{key}_std"] = float(np.std(values))
        summary[f"{key}_values"] = [float(value) for value in values]
    return summary


def _search_space(algorithm: str) -> dict[str, list]:
    spaces = {
        "random_forest": {"model__max_depth": [12, 30]},
        "extra_trees": {"model__max_depth": [12, 30]},
        "gradient_boosting": {"model__learning_rate": [0.03, 0.1]},
        "svm": {"model__C": [0.5, 2.0]},
        "kernel_ridge": {"model__alpha": [0.1, 1.0, 10.0]},
        "ridge": {"model__alpha": [0.1, 1.0, 10.0]},
        "elastic_net": {"model__alpha": [0.01, 0.1, 1.0]},
        "logistic_regression": {"model__C": [0.5, 2.0]},
        "mlp": {"model__alpha": [0.0001, 0.01]},
        "lightgbm": {"model__num_leaves": [15, 31]},
        "xgboost": {"model__max_depth": [4, 8]},
        "catboost": {"model__depth": [4, 8]},
    }
    return spaces.get(algorithm, {})


def _optimize(pipeline, smiles, targets, config):
    if not config.get("optimize_hyperparameters"):
        return pipeline, {"enabled": False}
    space = _search_space(config["algorithm"])
    if not space:
        return pipeline, {"enabled": True, "searched": False, "reason": "No safe comparison search space is defined."}
    folds = min(config["cross_validation_folds"], len(smiles))
    if config["task_type"] == "classification":
        _, counts = np.unique(targets, return_counts=True)
        folds = min(folds, int(counts.min()))
        splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=42)
        scoring = {"roc_auc": "roc_auc", "pr_auc": "average_precision", "f1": "f1_weighted", "accuracy": "accuracy"}.get(config.get("primary_metric"), "roc_auc")
    else:
        splitter = KFold(n_splits=folds, shuffle=True, random_state=42)
        scoring = {"rmse": "neg_root_mean_squared_error", "mae": "neg_mean_absolute_error", "r2": "r2"}.get(config.get("primary_metric"), "neg_root_mean_squared_error")
    search = GridSearchCV(pipeline, space, scoring=scoring, cv=splitter, refit=True, n_jobs=1, error_score="raise")
    search.fit(smiles, targets)
    return search.best_estimator_, {
        "enabled": True,
        "searched": True,
        "scoring": scoring,
        "candidate_count": len(search.cv_results_["params"]),
        "best_params": search.best_params_,
        "best_cv_score": float(search.best_score_),
    }


def _regression_metrics(actual, predicted):
    result = {
        "mae": float(mean_absolute_error(actual, predicted)),
        "rmse": float(math.sqrt(mean_squared_error(actual, predicted))),
    }
    if len(actual) >= 2:
        result["r2"] = float(r2_score(actual, predicted))
    else:
        result["r2"] = None
        result["r2_reason"] = "R²には2件以上のデータが必要です。"
    return result


def _classification_metrics(actual, predicted, probabilities, classes):
    average = "binary" if len(classes) == 2 else "weighted"
    positive = classes[1] if len(classes) == 2 else None
    result = {
        "accuracy": float(accuracy_score(actual, predicted)),
        "precision": float(precision_score(actual, predicted, average=average, pos_label=positive, zero_division=0)),
        "recall": float(recall_score(actual, predicted, average=average, pos_label=positive, zero_division=0)),
        "f1": float(f1_score(actual, predicted, average=average, pos_label=positive, zero_division=0)),
        "confusion_matrix": confusion_matrix(actual, predicted, labels=classes).tolist(),
        "classes": [str(item) for item in classes],
    }
    if len(classes) == 2 and len(np.unique(actual)) == 2:
        binary = (actual == positive).astype(int)
        result["roc_auc"] = float(roc_auc_score(binary, probabilities[:, 1]))
        result["pr_auc"] = float(average_precision_score(binary, probabilities[:, 1]))
    else:
        result["roc_auc"] = None
        result["pr_auc"] = None
        result["auc_reason"] = "ROC-AUCとPR-AUCは、Testに両方のクラスがある二値分類で表示します。"
    return result


def _nearest_similarities(test_smiles, train_smiles):
    train_fps = [ad_fingerprint(str(value)) for value in train_smiles]
    similarities = []
    for value in test_smiles:
        fingerprint = ad_fingerprint(str(value))
        scores = DataStructs.BulkTanimotoSimilarity(fingerprint, train_fps) if fingerprint else []
        similarities.append(max(scores) if scores else 0.0)
    return similarities


def train(config_path: Path, repository_root: Path):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    repository_root = repository_root.resolve()
    dataset_path = repository_root / config["dataset_path"]
    output_model = repository_root / config["output_model_path"]
    predictions_path = repository_root / config["predictions_path"]
    evaluation_path = repository_root / config["evaluation_path"]
    smiles, targets = _read_dataset(dataset_path, config["task_type"])
    train_smiles, test_smiles, train_targets, test_targets = _split(smiles, targets, config)
    print(f"QSAR split: training={len(train_smiles)} test={len(test_smiles)}", flush=True)

    pipeline = _pipeline(config)
    training_started = time.perf_counter()
    pipeline, optimization = _optimize(pipeline, train_smiles, train_targets, config)
    cv_metrics = _cross_validation(pipeline, train_smiles, train_targets, config)
    print("QSAR cross-validation completed", flush=True)
    pipeline.fit(train_smiles, train_targets)
    training_seconds = time.perf_counter() - training_started
    train_predictions = pipeline.predict(train_smiles)
    inference_started = time.perf_counter()
    test_predictions = pipeline.predict(test_smiles)
    inference_seconds = time.perf_counter() - inference_started
    similarities = _nearest_similarities(test_smiles, train_smiles)

    if config["task_type"] == "regression":
        training_metrics = _regression_metrics(train_targets, train_predictions)
        test_metrics = _regression_metrics(test_targets, test_predictions)
        probabilities = None
        residual_values = np.asarray(test_targets, dtype=float) - np.asarray(test_predictions, dtype=float)
        counts, edges = np.histogram(residual_values, bins=min(10, max(2, len(residual_values))))
        plots = {
            "actual_vs_predicted": [[float(actual), float(predicted)] for actual, predicted in zip(test_targets, test_predictions)],
            "residuals": [[float(predicted), float(actual - predicted)] for actual, predicted in zip(test_targets, test_predictions)],
            "residual_distribution": [[float((edges[index] + edges[index + 1]) / 2), int(count)] for index, count in enumerate(counts)],
        }
    else:
        train_probabilities = pipeline.predict_proba(train_smiles)
        probabilities = pipeline.predict_proba(test_smiles)
        classes = pipeline.named_steps["model"].classes_
        training_metrics = _classification_metrics(train_targets, train_predictions, train_probabilities, classes)
        test_metrics = _classification_metrics(test_targets, test_predictions, probabilities, classes)
        plots = {}
        if len(classes) == 2 and len(np.unique(test_targets)) == 2:
            binary = (test_targets == classes[1]).astype(int)
            false_positive, true_positive, _ = roc_curve(binary, probabilities[:, 1])
            precision, recall, _ = precision_recall_curve(binary, probabilities[:, 1])
            plots = {
                "roc_curve": [[float(x), float(y)] for x, y in zip(false_positive, true_positive)],
                "precision_recall_curve": [[float(x), float(y)] for x, y in zip(recall, precision)],
            }

    metrics = {
        "training": training_metrics,
        "test": test_metrics,
        "cross_validation": cv_metrics,
        "training_data_count": len(train_smiles),
        "test_data_count": len(test_smiles),
        "training_seconds": training_seconds,
        "inference_seconds": inference_seconds,
        "inference_molecules_per_second": len(test_smiles) / inference_seconds if inference_seconds else None,
        "applicability_domain_fraction": float(np.mean(np.asarray(similarities) >= config["ad_threshold"])),
        "performance_gap": {
            key: float(training_metrics[key] - test_metrics[key])
            for key in training_metrics.keys() & test_metrics.keys()
            if isinstance(training_metrics[key], (int, float)) and isinstance(test_metrics[key], (int, float))
        },
        "hyperparameter_optimization": optimization,
    }
    metadata = {
        **config,
        "library_versions": {
            "python": platform.python_version(),
            "rdkit": rdkit.__version__,
            "scikit_learn": sklearn.__version__,
            "joblib": joblib.__version__,
            "lightgbm": _package_version("lightgbm"),
            "xgboost": _package_version("xgboost"),
            "catboost": _package_version("catboost"),
        },
        "ad_method": "Morgan radius 2 / 2048-bit nearest-neighbour Tanimoto",
    }
    bundle = {
        "pipeline": pipeline,
        "metadata": metadata,
        "training_smiles": train_smiles.tolist(),
    }
    output_model.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, output_model)

    predictions_path.parent.mkdir(parents=True, exist_ok=True)
    with predictions_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        header = ["Split", "SMILES", "Actual", "Predicted", "NearestTanimoto", "In_AD"]
        if probabilities is not None:
            header.append("PositiveClassProbability")
        writer.writerow(header)
        for index, (value, actual, predicted, similarity) in enumerate(zip(test_smiles, test_targets, test_predictions, similarities)):
            row = ["test", value, actual, predicted, similarity, similarity >= config["ad_threshold"]]
            if probabilities is not None:
                row.append(float(probabilities[index, 1]) if probabilities.shape[1] == 2 else "")
            writer.writerow(row)
    evaluation_path.write_text(json.dumps({"metrics": metrics, "metadata": metadata, "plots": plots}, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("QSAR model and evaluation artifacts saved", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    arguments = parser.parse_args()
    train(arguments.config, arguments.repository_root)


if __name__ == "__main__":
    main()
