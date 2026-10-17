from pathlib import Path

import pytest

from app.services.results import ResultValidationError, validate_and_finalize_csv


def test_valid_csv_is_atomically_finalized(tmp_path: Path):
    source = tmp_path / "output" / "raw.csv"
    destination = tmp_path / "output" / "final.csv"
    source.parent.mkdir()
    source.write_text("SMILES,NLL\nCCO,1.2\n", encoding="utf-8")

    assert validate_and_finalize_csv(source, destination, tmp_path) == 1
    assert destination.read_text(encoding="utf-8") == "SMILES,NLL\nCCO,1.2\n"
    assert not source.exists()


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("", "empty"),
        ("NLL\n1.2\n", "SMILES"),
        ("SMILES,NLL\n", "no data rows"),
    ],
)
def test_invalid_csv_is_rejected(tmp_path: Path, content: str, message: str):
    source = tmp_path / "output" / "raw.csv"
    source.parent.mkdir()
    source.write_text(content, encoding="utf-8")

    with pytest.raises(ResultValidationError, match=message):
        validate_and_finalize_csv(
            source, tmp_path / "output" / "final.csv", tmp_path
        )


def test_symlink_and_path_escape_are_rejected(tmp_path: Path):
    outside = tmp_path.parent / f"{tmp_path.name}-outside.csv"
    outside.write_text("SMILES\nCCO\n", encoding="utf-8")
    source = tmp_path / "raw.csv"
    source.symlink_to(outside)

    with pytest.raises(ResultValidationError, match="escapes"):
        validate_and_finalize_csv(source, tmp_path / "final.csv", tmp_path)

    inside = tmp_path / "inside.csv"
    inside.write_text("SMILES\nCCO\n", encoding="utf-8")
    with pytest.raises(ResultValidationError, match="destination escapes"):
        validate_and_finalize_csv(inside, outside, tmp_path)

