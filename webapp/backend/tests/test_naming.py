from datetime import datetime, timezone
from uuid import UUID

import pytest

from app.services.naming import confined_path, result_filename, slugify_comment


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (" kinase A run1 ", "kinase_A_run1"),
        ("日本語 コメント", "日本語_コメント"),
        ("../secret", "secret"),
        ("/tmp/absolute", "tmp_absolute"),
        ("a\\b/c", "a_b_c"),
        ("\x00\x01", "result"),
        ("😀", "result"),
        ("...", "result"),
        ("", "result"),
        ("Ａ Ｂ", "A_B"),
    ],
)
def test_slugify_comment(value, expected):
    assert slugify_comment(value) == expected


def test_slug_is_bounded_and_ascii_mode_is_optional():
    assert len(slugify_comment("a" * 100)) == 50
    assert slugify_comment("日本語 abc", ascii_only=True) == "abc"


def test_result_filename_uses_configured_timezone_and_optional_collision_suffix():
    created = datetime(2026, 10, 16, 5, 30, tzinfo=timezone.utc)
    job_id = UUID("12345678-1234-5678-1234-567812345678")

    assert result_filename(
        comment_slug="kinase", created_at=created, timezone_name="Asia/Tokyo"
    ) == "202610161430_kinase.csv"
    assert result_filename(
        comment_slug="kinase",
        created_at=created,
        timezone_name="Asia/Tokyo",
        job_id=job_id,
    ) == "202610161430_kinase_12345678.csv"
    assert result_filename(
        comment_slug="comment", created_at=created, timezone_name="Asia/Tokyo",
        tab_type="MolGen", model_type="Mol2MolHigh"
    ) == "20261016_MolGen_Mol2MolHigh_comment.csv"


def test_confined_path_rejects_escape(tmp_path):
    assert confined_path(tmp_path, "jobs/one") == tmp_path / "jobs" / "one"
    with pytest.raises(ValueError):
        confined_path(tmp_path, "../outside")
