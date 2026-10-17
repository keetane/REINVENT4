# Phase 2 and Phase 3 implementation record

Date: 2026-10-16 (Asia/Tokyo)

## Phase 2

The API now persists jobs in SQLite and returns HTTP 202 without running computation in the request
process. Each UUID job receives isolated `input`, `config`, `logs`, and `output` directories.
Workers claim queued jobs with `BEGIN IMMEDIATE` plus a conditional status update, so two workers
cannot claim the same row.

The controlled fake subprocess covers the full lifecycle without model dependencies: startup,
bounded logs, output validation, atomic CSV finalization, nonzero exit, missing output, malformed
output, cancellation by process group, API disconnect, API restart, and stale-worker
reconciliation. A short reconciliation grace protects the interval between a claim and PID
registration; idle workers repeat reconciliation so genuinely abandoned jobs become
`interrupted` without automatic replay.

Comment slugs use Unicode NFKC, remove path/control characters, preserve safe Unicode letters and
digits, and are bounded to 50 characters. The final path is confined to the UUID job directory and
uses the configured `Asia/Tokyo` timezone.

## Phase 3

The REINVENT adapter uses a fixed model allowlist and validates device, sampling strategy, count,
temperature, and boolean flags. It writes a job-specific TOML snapshot and invokes:

```text
<configured-python> -m reinvent -l <job>/logs/reinvent.log <job>/config/sampling.toml
```

The subprocess receives an argument array with `shell=False`, a job-specific cwd, a minimal
environment, and an explicit server-configured source root via `PYTHONPATH`. User comments, model
paths, and arbitrary command arguments never enter the command line. Runtime/model/config/argv
metadata are saved in SQLite, while environment contents are not.

The existing `r4` editable mapping remains stale and was not modified. The adapter bypasses that
mapping explicitly and was verified with the current repository source, Python 3.10.17, REINVENT
4.5.11 metadata, PyTorch 2.6.0, and RDKit 2024.09.6.

## Verification

- Routine backend suite: 89 passed, 1 model-backed test skipped.
- Real LibInvent smoke: 1 passed in 6.38 seconds, two output molecules, exit code 0.
- Separate-process API restart scenario: job restored as `completed`, exit code 0, four-row CSV.
- Frontend Phase 1 suite remains unchanged: unit, production build, and Playwright E2E passed.
- npm audit: zero known vulnerabilities at the time of verification.

One deprecation warning remains in FastAPI/Starlette's `TestClient` compatibility layer; runtime
API and worker behavior are unaffected.

## Deferred work

- Phase 4 must validate and canonicalize parent structures with RDKit before accepting real jobs.
- Ketcher, history/progress UI, and polling recovery belong to Phase 4.
- Result query, descriptor normalization, filtering, statistics, and CSV/SDF downloads belong to
  Phase 5.
- Authentication is still absent, so services remain restricted to `127.0.0.1`.
