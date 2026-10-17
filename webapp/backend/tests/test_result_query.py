from pathlib import Path

from app.services.result_query import (
    ResultFilters,
    export_csv,
    export_sdf,
    load_result_dataset,
    query_rows,
    scatter_for,
    statistics_for,
)


def _csv(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "result.csv"
    path.write_text(content, encoding="utf-8")
    return path


def test_aliases_computed_descriptors_and_invalid_smiles(tmp_path: Path):
    dataset = load_result_dataset(
        _csv(
            tmp_path,
            "SMILES,CtRotBonds,RotB,NLL\nCCO,4,2,1.2\nnot-smiles,1,1,2.0\nCCCC,3,3,nan\n",
        )
    )

    assert dataset.total_rows == 3
    assert len(dataset.rows) == 2
    assert dataset.invalid_smiles == 1
    assert dataset.rows[0].values["RotB"] == 2
    assert dataset.rows[0].values["MW"] is not None
    assert 0 < dataset.rows[0].values["QED"] <= 1
    assert dataset.rows[0].values["SAScore"] is not None
    assert dataset.rows[0].values["HeavyAtomCount"] == 3
    assert dataset.rows[0].values["LipinskiViolations"] == 0
    assert dataset.rows[1].values["NLL"] is None
    assert "RotB aliases conflict" in dataset.warnings[0]


def test_filters_are_inclusive_and_sort_is_allowlisted(tmp_path: Path):
    dataset = load_result_dataset(_csv(tmp_path, "SMILES,MW\nCCO,40\nCCCC,60\nCCCCC,80\n"))
    rows = query_rows(
        dataset,
        ResultFilters(ranges={"MW": (40, 80)}, sort_by="MW", sort_direction="desc"),
    )

    assert [row.values["MW"] for row in rows] == [80, 60, 40]


def test_statistics_and_exports_share_the_same_filtered_rows(tmp_path: Path):
    dataset = load_result_dataset(
        _csv(tmp_path, "SMILES,MW,Name\nCCO,40,normal\nCCCC,60,=SUM(A1:A2)\nCCCCC,80,last\n")
    )
    rows = query_rows(dataset, ResultFilters(ranges={"MW": (50, 80)}))
    summaries, histograms = statistics_for(dataset, rows)
    csv_text = export_csv(dataset.columns, rows)
    sdf_text = export_sdf(rows)

    assert summaries["MW"]["count"] == 2
    assert summaries["MW"]["minimum"] == 60
    assert sum(item["count"] for item in histograms["MW"]) == 2
    assert csv_text.count("\n") == 3
    assert "'=SUM(A1:A2)" in csv_text
    assert sdf_text.count("$$$$") == 2

    points, sampled = scatter_for(dataset, rows, "MW", "LogP", limit=1)
    assert len(points) == 1
    assert points[0]["row_number"] == 2
    assert sampled is True


def test_parent_and_scaffold_similarity_use_the_job_parent(tmp_path: Path):
    dataset = load_result_dataset(
        _csv(tmp_path, "SMILES\nc1ccccc1\nC1CCCCC1\n"),
        parent_smiles="c1ccccc1",
    )

    assert dataset.rows[0].values["ParentSimilarity"] == 1.0
    assert dataset.rows[0].values["ScaffoldSimilarity"] == 1.0
    assert dataset.rows[1].values["ParentSimilarity"] < 1.0
    assert dataset.rows[1].values["ScaffoldSimilarity"] < 1.0
