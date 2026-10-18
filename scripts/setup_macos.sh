#!/bin/zsh
set -euo pipefail

ROOT_DIR="${0:A:h:h}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

cd "$ROOT_DIR"
mkdir -p .runtimes

"$PYTHON_BIN" -m venv .runtimes/reinvent
source .runtimes/reinvent/bin/activate
python -m pip install --upgrade pip
python install.py mac -d none -e
python -m pip install 'lightgbm>=4,<5' 'xgboost>=2,<4' 'catboost>=1.2,<2'
deactivate

"$PYTHON_BIN" -m venv webapp/backend/.venv
webapp/backend/.venv/bin/python -m pip install --upgrade pip
webapp/backend/.venv/bin/pip install -r webapp/backend/requirements.lock
webapp/backend/.venv/bin/pip install -e webapp/backend --no-deps

"$PYTHON_BIN" -m venv .runtimes/chemprop-v2
.runtimes/chemprop-v2/bin/python -m pip install --upgrade pip
.runtimes/chemprop-v2/bin/python -m pip install 'chemprop==2.3.1'

"$PYTHON_BIN" -m venv .runtimes/qsartuna
.runtimes/qsartuna/bin/python -m pip install --upgrade pip
.runtimes/qsartuna/bin/python -m pip install 'git+https://github.com/MolecularAI/QSARtuna.git@664f73e3a9f68aaf788ee6a92a4889778ddeaf84' 'optuna>=3.1,<4' 'sqlalchemy==1.4.46' 'mapchiral==0.0.5' 'peptides==0.5.0'

cd webapp/frontend
npm ci
cd "$ROOT_DIR"
python scripts/check_environment.py
