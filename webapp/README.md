# REINVENT4 Web application

This directory is the local-first replacement for the legacy Streamlit UI. It has three
independently started processes:

- FastAPI serves the HTTP API on `127.0.0.1:8000`.
- The worker is a separate process and atomically claims queued jobs.
- React/Vite serves the browser UI on `127.0.0.1:5173`.

The API and worker both apply the same numbered SQLite migrations. The default database is
`<repository>/jobs.db`, job artifacts will live under `<repository>/jobs`, and both paths can be
overridden with environment variables documented in `backend/.env.example`.

## Backend

```bash
cd webapp/backend
python -m venv .venv
.venv/bin/pip install -r requirements.lock
.venv/bin/pip install -e . --no-deps
.venv/bin/uvicorn app.main:app --host 127.0.0.1 --port 8000
```

In a second terminal, start the independent worker:

```bash
cd webapp/backend
.venv/bin/python -m worker
```

For a non-blocking initialization check, use `.venv/bin/python -m worker --once`.

The default test runner deliberately uses a controlled fake subprocess instead of REINVENT. It exercises process
separation, cancellation, logging, output validation, and crash recovery without depending on a
model runtime. It remains the default runner for routine tests.

The API endpoints are:

- `POST /api/parents/validate`
- `POST /api/jobs`
- `GET /api/jobs`
- `GET /api/jobs/{job_id}`
- `POST /api/jobs/{job_id}/cancel`
- `GET /api/jobs/{job_id}/logs?tail=100`
- `GET /api/jobs/{job_id}/results/schema`
- `GET /api/jobs/{job_id}/results/molecules`
- `GET /api/jobs/{job_id}/results/statistics`
- `GET /api/jobs/{job_id}/results/download.csv`
- `GET /api/jobs/{job_id}/results/download.sdf`
- `POST /api/qsar-data/inspect`
- `POST /api/qsar-jobs`
- `POST /api/qsar-comparison-jobs`
- `GET /api/qsar-comparisons/{job_id}`
- `POST /api/qsar-comparisons/{job_id}/register`
- `GET/DELETE /api/qsar-models[/{model_id}]`
- `GET /api/qsar-models/{model_id}/artifacts/{model|dataset|predictions|evaluation}`
- `POST /api/qsar/predict`
- `POST /api/rl-jobs`
- `POST /api/rl-jobs/{job_id}/progress`

## QSAR and Reinforcement Learning

The QSAR tab accepts CSV data, validates and canonicalizes SMILES with RDKit, records rejected
rows, resolves duplicate conflicts according to the selected policy, and optionally converts M,
mM, μM/uM, nM, or pM concentrations to pActivity. Random Forest, SVM, Gradient Boosting, Extra
Trees, Kernel Ridge, Gaussian Process, Ridge, Elastic Net, Logistic Regression, MLP, LightGBM,
XGBoost, and CatBoost are available for their supported regression/classification tasks.
Fingerprints, RDKit descriptors, train-only preprocessing, random/scaffold split,
cross-validation, test metrics, timing, and Morgan-neighbour applicability-domain results are
stored with the trusted model artifacts.

Model Comparison applies the same deterministic split to every selected algorithm. The recommended
model is selected only from training cross-validation results; the held-out Test set is reported
afterward and is not used for model selection. Successful candidates can be registered separately.

Chemprop v2.3.1 and QSARtuna 4.0.1 run in separate virtual environments under `.runtimes/` so their
PyTorch, RDKit, and optimization dependencies do not replace the REINVENT runtime. Chemprop uses
its official v2 command-line training and prediction interfaces. QSARtuna uses its official Python
optimization/build API. Their saved models are connected to REINVENT staged learning through the
official `ExternalProcess` scoring component.

The RL tab creates an independent `staged_learning` job. It can start from a standard Prior, a TL
model, or an earlier RL checkpoint. The server snapshots all selected generation/QSAR model IDs and
paths when the job is queued. A repository scoring plugin named `QsarJoblib` loads only trusted
server-created model bundles, performs batched predictions, and applies the stored Morgan nearest-
neighbour AD threshold. Invalid or non-finite predictions are never converted into a high score.
REINVENT's official aggregation, transforms, DAP learning strategy, responder, diversity filter,
summary CSV, and checkpoint handling remain in use.

With `REINVENT_WEB_WORKER_RUNNER=reinvent`, QSAR training and prediction use the configured REINVENT
runtime because it already contains RDKit 2024.09.6 and scikit-learn 1.6.1. The API virtualenv does
not deserialize QSAR models. Model downloads are trusted artifacts and must never be uploaded back
and loaded from an untrusted source.

## REINVENT adapter

The adapter does not repair or depend on the stale editable installation in the existing
`r4` environment. It invokes the configured Python as a module, sets `PYTHONPATH` to the configured
live source tree, uses the job directory as cwd, and only accepts server-defined model keys.

For the currently verified local runtime:

```bash
cd webapp/backend
export REINVENT_WEB_WORKER_RUNNER=reinvent
export REINVENT_RUNTIME_PYTHON=/Users/keetane/miniforge3/envs/r4/bin/python
export REINVENT_SOURCE_ROOT=/Users/keetane/Documents/modeling/REINVENT4
.venv/bin/python -m worker
```

The currently allowed model keys are `mol2mol_high`, `mol2mol_medium`, `mol2mol_mmp`,
`mol2mol_scaffold`, `mol2mol_scaffold_generic`, `mol2mol_similarity`, `reinvent`, `libinvent`, and
`linkinvent`. A client supplies a key, never an executable or model path. A minimal job body is:

```json
{
  "parent_smiles": "N#CC[C@H](C1CCCC1)[*:1]",
  "comment": "libinvent smoke",
  "run_conditions": {
    "model_key": "libinvent",
    "num_smiles": 2,
    "device": "cpu"
  }
}
```

The real model smoke is intentionally separate from the routine suite:

```bash
RUN_REINVENT_SMOKE=1 \
REINVENT_SMOKE_SOURCE=/Users/keetane/Documents/modeling/REINVENT4 \
REINVENT_SMOKE_PYTHON=/Users/keetane/miniforge3/envs/r4/bin/python \
.venv/bin/pytest -m reinvent -vv
```

An RL smoke should use one step and a small batch, then verify all three artifacts: the checkpoint,
`rl-summary_1.csv`, and a sampling CSV created from the checkpoint. The test is intentionally not in
the routine suite because the Prior is roughly 22 MB and the run imports PyTorch and all scoring
plugins.

No stable per-step progress signal was observed in the sampling logs, so running jobs expose phase
messages only. A percentage is stored only after validated completion.

## Frontend

Use Node 24.14.1 or newer; Ketcher 3.18 declares that engine range. Ketcher is loaded only when
the structure editor is opened, runs locally, and the text SMILES path remains available if it
cannot initialize.

```bash
cd webapp/frontend
npm ci
npm run dev
```

Open `http://127.0.0.1:5173`. The Vite development proxy forwards `/api` to the loopback API.
The language selector stores Japanese/English preference in local storage. Navigation and Tutorial
content are translated from `src/i18n/`; feature-form translation is being migrated incrementally.

## Clean setup

The supported setup scripts install REINVENT, the web backend, Chemprop v2, QSARtuna, and frontend
dependencies. They deliberately do not install OpenEye or ROCS.

```bash
# Apple Silicon macOS
bash scripts/setup_macos.sh

# Ubuntu (CPU by default; set REINVENT_PROCESSOR_TYPE before running for CUDA)
bash scripts/setup_ubuntu.sh

python scripts/check_environment.py
```

Every required runtime must report `[OK]`. The checker exits non-zero when a runtime or command is
missing. See the in-app Tutorial for the same instructions in Japanese and English.

## Result semantics

`MW`, `RotB`, `HBD`, `HBA`, `LogP`, and `TPSA` can be filtered with inclusive bounds and combined
with AND. Missing/non-finite values are excluded only when their descriptor is filtered. When a
required descriptor is absent, the server computes it with RDKit 2024.09.6 using `MolWt`,
`NumRotatableBonds`, `NumHDonors`, `NumHAcceptors`, `MolLogP`, and `CalcTPSA`. Existing values are
preserved. `RotB` uses the precedence `RotB`, `RoB`, `CtRotB`, `CtRotBonds`,
`Rotatable_Bonds`, and reports conflicting aliases.

Statistics and CSV/SDF downloads use the same server-side filter service as the paginated browser
view. Invalid SMILES are counted and skipped without hiding valid rows. CSV cells that begin with a
spreadsheet formula marker are neutralized unless they are valid numbers. Tanimoto values are never
invented; an existing `Tanimoto_Similarity_to_Parent` column is used as supplied.

## Verification

```bash
cd webapp/backend
.venv/bin/pytest

cd ../frontend
npm test
npm run build
npm run e2e
```

No service in this phase binds to `0.0.0.0`. Do not expose the unauthenticated application to a
LAN or the internet.

See [operations.md](docs/operations.md) for startup, shutdown, recovery, backup, staged migration,
and rollback. See [security-review.md](docs/security-review.md) and
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) before changing the deployment boundary. See
[qsar-rl.md](docs/qsar-rl.md) for model semantics, supported versions, and operating limits.
