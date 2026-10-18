from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

from worker.qsar_training import train


def _ranking_value(result: dict, metric: str) -> float:
    value = result.get("metrics", {}).get("cross_validation", {}).get(metric)
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        return math.inf if metric in {"rmse", "mae"} else -math.inf
    return float(value)


def compare(config_path: Path, repository_root: Path) -> None:
    request = json.loads(config_path.read_text(encoding="utf-8"))
    results: list[dict] = []
    total = len(request["models"])
    for index, item in enumerate(request["models"], start=1):
        algorithm = item["algorithm"]
        print(f"QSAR comparison {index}/{total}: {algorithm}", flush=True)
        started = time.perf_counter()
        try:
            train(repository_root / item["config_path"], repository_root)
            evaluation = json.loads((repository_root / item["evaluation_path"]).read_text(encoding="utf-8"))
            results.append({
                "algorithm": algorithm,
                "status": "completed",
                "duration_seconds": time.perf_counter() - started,
                "metrics": evaluation["metrics"],
                "plots": evaluation.get("plots", {}),
                "model_path": item["output_model_path"],
                "predictions_path": item["predictions_path"],
                "evaluation_path": item["evaluation_path"],
            })
        except Exception as error:
            results.append({
                "algorithm": algorithm,
                "status": "failed",
                "duration_seconds": time.perf_counter() - started,
                "error": f"{type(error).__name__}: {error}",
            })
            print(f"QSAR comparison failed for {algorithm}: {type(error).__name__}: {error}", flush=True)

    successful = [result for result in results if result["status"] == "completed"]
    if not successful:
        raise RuntimeError("all algorithms in the comparison failed")
    metric = request["primary_metric"]
    reverse = metric not in {"rmse", "mae"}
    recommended = sorted(successful, key=lambda result: _ranking_value(result, metric), reverse=reverse)[0]
    summary = {
        "primary_metric": metric,
        "recommended_algorithm": recommended["algorithm"],
        "selection_basis": "training-only cross-validation",
        "results": results,
        "successful_count": len(successful),
        "failed_count": len(results) - len(successful),
    }
    output = repository_root / request["summary_path"]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Recommended algorithm: {recommended['algorithm']} ({metric}, cross-validation)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    arguments = parser.parse_args()
    compare(arguments.config, arguments.repository_root.resolve())


if __name__ == "__main__":
    main()
