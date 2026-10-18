from pathlib import Path

from app.config import Settings


def test_runtime_python_paths_keep_virtual_environment_symlinks(monkeypatch, tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    runtime = repository / ".runtimes" / "chemprop-v2" / "bin" / "python"
    runtime.parent.mkdir(parents=True)
    target = tmp_path / "shared" / "python"
    target.parent.mkdir()
    target.touch()
    runtime.symlink_to(target)

    monkeypatch.setenv("REINVENT_WEB_REPOSITORY_ROOT", str(repository))
    monkeypatch.setenv("CHEMPROP_V2_PYTHON", str(runtime))

    settings = Settings.from_env()

    assert settings.chemprop_runtime_python == runtime
    assert settings.chemprop_runtime_python != runtime.resolve()
