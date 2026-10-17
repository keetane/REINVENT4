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
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) before changing the deployment boundary.
