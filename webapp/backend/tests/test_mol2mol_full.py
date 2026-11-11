from __future__ import annotations

import csv
from pathlib import Path

from worker.mol2mol_full import SOURCE_MODEL_COLUMN, combine_results


def _output(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SMILES", "NLL"])
        writer.writeheader()
        writer.writerows(rows)


def test_combines_results_with_source_model_as_last_column(tmp_path: Path):
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    destination = tmp_path / "combined.csv"
    _output(first, [{"SMILES": "CCO", "NLL": "1.0"}])
    _output(second, [{"SMILES": "CCC", "NLL": "2.0"}])
    runs = [
        {"model_key": "mol2mol_high", "raw_output": str(first)},
        {"model_key": "mol2mol_medium", "raw_output": str(second)},
    ]

    count = combine_results(runs, destination, tmp_path)

    with destination.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    assert count == 2
    assert reader.fieldnames == ["SMILES", "NLL", SOURCE_MODEL_COLUMN]
    assert [row[SOURCE_MODEL_COLUMN] for row in rows] == [
        "mol2mol_high",
        "mol2mol_medium",
    ]
