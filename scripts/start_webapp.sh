#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BACKEND_DIR="$ROOT_DIR/webapp/backend"
FRONTEND_DIR="$ROOT_DIR/webapp/frontend"
RUN_DIR="$ROOT_DIR/webapp/.run"
BACKEND_PYTHON="$BACKEND_DIR/.venv/bin/python"
REINVENT_PYTHON="${REINVENT_RUNTIME_PYTHON:-$ROOT_DIR/.runtimes/reinvent/bin/python}"

if [[ ! -x "$REINVENT_PYTHON" && -x "$HOME/miniforge3/envs/r4/bin/python" ]]; then
  REINVENT_PYTHON="$HOME/miniforge3/envs/r4/bin/python"
fi

[[ -x "$BACKEND_PYTHON" ]] || { echo "Missing backend runtime. Run the setup script first." >&2; exit 1; }
[[ -x "$REINVENT_PYTHON" ]] || { echo "Missing REINVENT runtime. Run the setup script first." >&2; exit 1; }
command -v npm >/dev/null || { echo "npm is required." >&2; exit 1; }
command -v curl >/dev/null || { echo "curl is required." >&2; exit 1; }

mkdir -p "$RUN_DIR"

running() {
  [[ -f "$1" ]] || return 1
  local pid
  pid="$(<"$1")"
  [[ "$pid" =~ ^[0-9]+$ ]] && kill -0 "$pid" 2>/dev/null
}

if ! curl -fsS --max-time 1 http://127.0.0.1:8000/api/health >/dev/null 2>&1 && ! running "$RUN_DIR/api.pid"; then
  (cd "$BACKEND_DIR" && nohup env PYTHONPATH="$BACKEND_DIR" REINVENT_WEB_REPOSITORY_ROOT="$ROOT_DIR" REINVENT_WEB_DATABASE_PATH="$ROOT_DIR/jobs.db" REINVENT_WEB_JOBS_DIR="$ROOT_DIR/jobs" REINVENT_WEB_WORKER_RUNNER=reinvent REINVENT_RUNTIME_PYTHON="$REINVENT_PYTHON" REINVENT_SOURCE_ROOT="$ROOT_DIR" "$BACKEND_PYTHON" -m uvicorn app.main:app --host 127.0.0.1 --port 8000 >>"$RUN_DIR/api.log" 2>&1 & echo $! >"$RUN_DIR/api.pid")
fi

if ! running "$RUN_DIR/worker.pid"; then
  (cd "$BACKEND_DIR" && nohup env PYTHONPATH="$BACKEND_DIR" REINVENT_WEB_REPOSITORY_ROOT="$ROOT_DIR" REINVENT_WEB_DATABASE_PATH="$ROOT_DIR/jobs.db" REINVENT_WEB_JOBS_DIR="$ROOT_DIR/jobs" REINVENT_WEB_WORKER_RUNNER=reinvent REINVENT_RUNTIME_PYTHON="$REINVENT_PYTHON" REINVENT_SOURCE_ROOT="$ROOT_DIR" "$BACKEND_PYTHON" -m worker >>"$RUN_DIR/worker.log" 2>&1 & echo $! >"$RUN_DIR/worker.pid")
fi

if ! curl -fsS --max-time 1 http://127.0.0.1:5173/ >/dev/null 2>&1 && ! running "$RUN_DIR/frontend.pid"; then
  (cd "$FRONTEND_DIR" && nohup npm run dev -- --host 127.0.0.1 --port 5173 >>"$RUN_DIR/frontend.log" 2>&1 & echo $! >"$RUN_DIR/frontend.pid")
fi

for _ in $(seq 1 120); do
  if curl -fsS --max-time 1 http://127.0.0.1:8000/api/health >/dev/null 2>&1 && curl -fsS --max-time 1 http://127.0.0.1:5173/ >/dev/null 2>&1; then
    echo "Molecular Generator is ready: http://127.0.0.1:5173/"
    if [[ "${OPEN_BROWSER:-1}" == "1" ]]; then
      if command -v open >/dev/null; then open http://127.0.0.1:5173/; elif command -v xdg-open >/dev/null; then xdg-open http://127.0.0.1:5173/ >/dev/null 2>&1 || true; fi
    fi
    exit 0
  fi
  sleep 1
done

echo "Startup timed out. See $RUN_DIR/api.log, worker.log, and frontend.log." >&2
exit 1
