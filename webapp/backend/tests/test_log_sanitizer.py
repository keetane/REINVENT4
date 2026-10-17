from app.services.log_sanitizer import sanitize_log_lines


def test_log_lines_hide_known_paths_and_secret_assignments(initialized_settings):
    settings = initialized_settings
    lines = [
        f"reading {settings.jobs_dir}/job/output.csv",
        "API_KEY=should-not-appear password:also-secret",
    ]

    sanitized = sanitize_log_lines(lines, settings)

    assert str(settings.jobs_dir) not in sanitized[0]
    assert "<jobs>" in sanitized[0]
    assert "should-not-appear" not in sanitized[1]
    assert "also-secret" not in sanitized[1]
