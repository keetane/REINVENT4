from worker.__main__ import initialize_worker


def test_worker_can_initialize_independently(settings):
    initialize_worker(settings)

    assert settings.database_path.is_file()
    assert settings.jobs_dir.is_dir()

