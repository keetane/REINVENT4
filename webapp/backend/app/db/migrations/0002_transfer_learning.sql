ALTER TABLE jobs ADD COLUMN job_type TEXT NOT NULL DEFAULT 'generation'
    CHECK (job_type IN ('generation', 'transfer_learning'));
ALTER TABLE jobs ADD COLUMN model_id TEXT;
ALTER TABLE jobs ADD COLUMN metrics_json TEXT;

CREATE INDEX jobs_model_id_status_idx ON jobs(model_id, status);

CREATE TABLE trained_models (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    model_type TEXT NOT NULL CHECK (model_type IN ('reinvent', 'mol2mol')),
    created_at TEXT NOT NULL,
    training_data_count INTEGER NOT NULL CHECK (training_data_count >= 1),
    validation_data_count INTEGER NOT NULL CHECK (validation_data_count >= 0),
    training_conditions_json TEXT NOT NULL,
    file_path TEXT NOT NULL UNIQUE,
    config_path TEXT NOT NULL,
    log_path TEXT NOT NULL,
    source_job_id TEXT NOT NULL UNIQUE,
    base_model_key TEXT NOT NULL,
    FOREIGN KEY (source_job_id) REFERENCES jobs(id)
);

CREATE INDEX trained_models_created_at_idx ON trained_models(created_at DESC);
