PRAGMA defer_foreign_keys = ON;

CREATE TABLE jobs_v4 (
    id TEXT PRIMARY KEY,
    status TEXT NOT NULL CHECK (status IN ('queued', 'running', 'completed', 'failed', 'cancelling', 'cancelled', 'interrupted')),
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
    progress_percent REAL CHECK (progress_percent IS NULL OR (progress_percent >= 0 AND progress_percent <= 100)),
    progress_message TEXT,
    job_dir TEXT NOT NULL UNIQUE,
    log_path TEXT,
    result_csv_path TEXT,
    result_sdf_path TEXT,
    error_summary TEXT,
    error_kind TEXT,
    cancel_requested_at TEXT,
    job_type TEXT NOT NULL DEFAULT 'generation' CHECK (
        job_type IN ('generation', 'transfer_learning', 'qsar_training', 'qsar_comparison', 'reinforcement_learning')
    ),
    model_id TEXT,
    metrics_json TEXT
);

INSERT INTO jobs_v4 SELECT * FROM jobs;
DROP TABLE jobs;
ALTER TABLE jobs_v4 RENAME TO jobs;
CREATE INDEX jobs_status_queued_at_idx ON jobs(status, queued_at);
CREATE INDEX jobs_created_at_idx ON jobs(created_at DESC);
CREATE INDEX jobs_model_id_status_idx ON jobs(model_id, status);

CREATE TABLE qsar_models_v4 (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    algorithm TEXT NOT NULL CHECK (
        algorithm IN (
            'random_forest', 'svm', 'gradient_boosting', 'extra_trees',
            'lightgbm', 'xgboost', 'catboost', 'kernel_ridge',
            'gaussian_process', 'ridge', 'elastic_net',
            'logistic_regression', 'mlp', 'chemprop_v2', 'qsartuna'
        )
    ),
    task_type TEXT NOT NULL CHECK (task_type IN ('regression', 'classification')),
    target_name TEXT NOT NULL,
    target_unit TEXT NOT NULL,
    created_at TEXT NOT NULL,
    training_data_count INTEGER NOT NULL CHECK (training_data_count >= 1),
    test_data_count INTEGER NOT NULL CHECK (test_data_count >= 1),
    metrics_json TEXT NOT NULL,
    metadata_json TEXT NOT NULL,
    file_path TEXT NOT NULL UNIQUE,
    dataset_path TEXT NOT NULL,
    predictions_path TEXT NOT NULL,
    evaluation_path TEXT NOT NULL,
    source_job_id TEXT NOT NULL,
    ad_threshold REAL NOT NULL CHECK (ad_threshold >= 0 AND ad_threshold <= 1),
    FOREIGN KEY (source_job_id) REFERENCES jobs(id)
);

INSERT INTO qsar_models_v4 SELECT * FROM qsar_models;
DROP TABLE qsar_models;
ALTER TABLE qsar_models_v4 RENAME TO qsar_models;
CREATE INDEX qsar_models_created_at_idx ON qsar_models(created_at DESC);
CREATE INDEX qsar_models_source_job_id_idx ON qsar_models(source_job_id);
