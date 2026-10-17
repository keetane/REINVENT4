# Operations guide

## Supported boundary

This release is a single-user, loopback-only application. Run the API on `127.0.0.1:8000` and the
frontend on `127.0.0.1:5173`. Do not bind either service to `0.0.0.0` or expose it to a LAN or the
internet. Authentication, authorization, TLS, CSRF policy, rate limiting, audit logging, retention,
and multi-user ownership are not implemented.

## Initial setup and startup

Create the backend environment and install the locked dependencies as shown in the root README.
Copy `backend/.env.example` into the process manager or shell environment and replace every
absolute-path placeholder. Keep `REINVENT_WEB_WORKER_RUNNER=fake` until the local API/worker
lifecycle is verified. For real jobs, set it to `reinvent` and configure the runtime Python and
live REINVENT4 source root explicitly.

Start, in order, one API process, one worker process, and the Vite frontend. Multiple workers are
safe for claiming, but the initial operating model uses one worker to avoid accidental contention
for CPU/GPU and model files. Confirm `GET http://127.0.0.1:8000/api/health` before creating jobs.

## Shutdown and restart

- Stop the frontend at any time; jobs are unaffected.
- Stop and restart the API at any time; the independent worker and child model process continue.
- Prefer stopping the worker while idle. A forced worker death leaves its job `running`; after the
  reconciliation grace period, the next worker marks an abandoned job `interrupted`. Jobs are not
  automatically replayed.
- To stop an active job deliberately, use the UI cancel action and wait for `cancelled`; it targets
  only that job's process group.

## Recovery and troubleshooting

1. Check the health endpoint, then confirm the worker process is running.
2. Inspect the job status and bounded log tail in the UI. API output redacts configured local paths
   and common secret assignments; the job's on-disk log remains an operator-only artifact.
3. A `failed` job retains the original exit code and a bounded diagnostic summary. A missing CSV is
   treated as output validation failure, not as the root cause by itself.
4. An `interrupted` job is intentionally not resumed. Create a new job after checking model inputs,
   runtime paths, disk space, and the original log.
5. If SQLite reports contention, ensure the database is on a local filesystem and only expected API
   and worker processes are using it. WAL and a busy timeout are enabled.

## Backup and restore

Stop job creation and wait for active work to finish. Stop the API and worker, then copy `jobs.db`
plus the entire `jobs/` tree as one backup set. Restore both to the same relative state while all
processes are stopped, verify ownership and free space, start the API so migrations run, then start
the worker. Never restore only the database or only the artifacts and assume they are consistent.

No automatic retention or deletion policy exists. Backups may contain molecular structures,
comments, logs, and model outputs and must receive the same access controls as the working copy.

## Staged migration and rollback

Run the new stack beside the existing Streamlit application on different loopback ports. Start with
fake jobs, then one approved two-molecule real smoke, then compare representative legacy results and
descriptor/filter behavior. Keep Streamlit unchanged and available until users accept the new flow.

Rollback is operational: stop the React, API, and worker processes and resume the existing
Streamlit command. The new implementation does not modify legacy application code or legacy result
files. Preserve `jobs.db` and `jobs/` for diagnosis; removing them is not part of rollback. If a
schema rollback is ever required, restore a matching database-and-artifact backup rather than
editing migration history.

## Deployment changes requiring a separate review

Any non-loopback exposure, multiple-user access, new upload support, a reverse proxy, background
service installation, credential handling, automated retention, or shared GPU scheduling changes
the threat model and requires explicit design and security review before rollout.
