from __future__ import annotations

import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True, type=Path)
    parser.add_argument("--evaluation", required=True, type=Path)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--training-count", required=True, type=int)
    parser.add_argument("--config", required=True, type=Path)
    arguments = parser.parse_args()
    arguments.model.write_bytes(b"trusted fake QSAR model")
    test_count = max(1, arguments.training_count // 5)
    train_count = arguments.training_count - test_count
    result = {
        "metrics": {
            "training": {"r2": 0.9, "mae": 0.2, "rmse": 0.3},
            "test": {"r2": 0.7, "mae": 0.4, "rmse": 0.5},
            "cross_validation": {"available": True, "folds": 3, "r2": 0.6},
            "training_data_count": train_count,
            "test_data_count": test_count,
        },
        "metadata": {**json.loads(arguments.config.read_text(encoding="utf-8")), "fake": True},
    }
    arguments.evaluation.write_text(json.dumps(result), encoding="utf-8")
    arguments.predictions.write_text(
        "Split,SMILES,Actual,Predicted,NearestTanimoto,In_AD\n"
        "test,CCO,5.0,5.1,0.8,True\n",
        encoding="utf-8",
    )
    print("fake QSAR training completed", flush=True)


if __name__ == "__main__":
    main()
