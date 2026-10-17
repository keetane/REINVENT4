CREATE TABLE jobs (
    id TEXT PRIMARY KEY,
    status TEXT NOT NULL CHECK (
        status IN (
            'queued', 'running', 'completed', 'failed',
            'cancelling', 'cancelled', 'interrupted'
        )
    ),
    created_at TEXT NOT NULL,
    queued_at TEXT NOT NULL,
    started_at TEXT,
    finished_at TEXT,
    updated_at TEXT NOT NULL,
    comment_original TEXT NOT NULL DEFAULT '',
    comment_slug TEXT NOT NULL,
    parent_smiles_input TEXT NOT NULL,
    parent_smiles_canonical TEXT,
    run_conditions_json TEXT NOT NULL,
    config_snapshot_path TEXT,
    model_path TEXT,
    command_metadata_json TEXT,
    worker_id TEXT,
    process_id INTEGER,
    exit_code INTEGER,
    progress_current INTEGER,
    progress_total INTEGER,
    progress_percent REAL CHECK (
        progress_percent IS NULL OR
        (progress_percent >= 0 AND progress_percent <= 100)
    ),
    progress_message TEXT,
    job_dir TEXT NOT NULL UNIQUE,
    log_path TEXT,
    result_csv_path TEXT,
    result_sdf_path TEXT,
    error_summary TEXT,
    error_kind TEXT,
    cancel_requested_at TEXT
);

CREATE INDEX jobs_status_queued_at_idx ON jobs(status, queued_at);
CREATE INDEX jobs_created_at_idx ON jobs(created_at DESC);

