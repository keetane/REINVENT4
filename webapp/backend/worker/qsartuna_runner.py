from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/qsartuna-matplotlib")

from apischema import serialize
import dill
from qsartuna.config import ModelMode, OptimizationDirection
from qsartuna.config.optconfig import (
    LogisticRegression,
    OptimizationConfig,
    RandomForestClassifier,
    RandomForestRegressor,
    Ridge,
    SVC,
    SVR,
)
from qsartuna.datareader import Dataset
from qsartuna.descriptors import ECFP, MACCS_keys
from qsartuna.three_step_opt_build_merge import build_best, buildconfig_best, optimize
from qsartuna.utils.preprocessing.splitter import Random, ScaffoldSplit


def main() -> None:
    # QSARtuna's model writer uses dill. Referencing imported enums/classes avoids
    # Python 3.12 attempting to recreate Enum internals when the model is loaded.
    dill.settings["byref"] = True
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--build-config", required=True)
    parser.add_argument("--task-type", choices=["regression", "classification"], required=True)
    parser.add_argument("--split-method", choices=["scaffold", "random"], required=True)
    parser.add_argument("--folds", type=int, required=True)
    parser.add_argument("--trials", type=int, required=True)
    parser.add_argument("--startup-trials", type=int, required=True)
    args = parser.parse_args()
    classification = args.task_type == "classification"
    algorithms = (
        [RandomForestClassifier.new(), SVC.new(), LogisticRegression.new()]
        if classification else [RandomForestRegressor.new(), SVR.new(), Ridge.new()]
    )
    splitter = ScaffoldSplit(random_state=42) if args.split_method == "scaffold" else Random(fraction=0.2, seed=42)
    config = OptimizationConfig(
        data=Dataset(training_dataset_file=args.train, test_dataset_file=args.test, input_column="Canonical_SMILES", response_column="Target"),
        descriptors=[ECFP.new(), MACCS_keys.new()],
        algorithms=algorithms,
        settings=OptimizationConfig.Settings(
            mode=ModelMode.CLASSIFICATION if classification else ModelMode.REGRESSION,
            optimization_split_strategy=splitter,
            direction=OptimizationDirection.MAXIMIZATION,
            n_trials=args.trials,
            n_startup_trials=min(args.startup_trials, args.trials),
            n_splits=args.folds,
            n_replicates=1,
            n_jobs=1,
            random_seed=42,
        ),
    )
    study = optimize(config, study_name="reinvent4-web-qsartuna")
    build_config = buildconfig_best(study)
    Path(args.build_config).write_text(json.dumps(serialize(build_config), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    build_best(build_config, args.model)


if __name__ == "__main__":
    main()
