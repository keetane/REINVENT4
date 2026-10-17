# Phase 0 investigation

Date: 2026-10-16 (Asia/Tokyo)

## Scope and repository state

The live repository at `/Users/keetane/Documents/modeling/REINVENT4` was inspected read-only before
application changes. `main` matched `origin/main` at `f12dbb3` and the tracked worktree was clean.
No repository-specific `AGENTS.md` was present between the repository and filesystem root.

## Existing Streamlit behavior

`reinventer.py` recursively lists CSV files under the current results directory, reads UTF-8 CSV,
sorts and filters known descriptors, renders the first 100 structures, and exports the current
filtered rows as CSV or SDF. It assumes the selected sort column and `SMILES` exist, so a missing
column, malformed CSV, invalid molecule, or disappearing file can currently fail the page.

`pages/auto_R.py` accepts a parent SMILES (or resolves a compound name through PubChem), performs
RECAP decomposition, and runs Mol2Mol, LibInvent, or LinkInvent sampling synchronously in the
Streamlit request process. It then calculates RDKit descriptors and rewrites the raw CSV in place.

## Observed REINVENT command and environment

The code constructs an argument array and calls:

```text
reinvent -l <working-directory>/results/log/sampling.log \
  <working-directory>/toml/sampling.toml
```

`shell=True` is not used. The subprocess inherits Streamlit's cwd and full environment. The
documented launch is from the repository root (`streamlit run reinventer.py`), so the observed cwd
contract is the repository root. Generated TOML selects `run_type = "sampling"` and supplies
`device`, model, output CSV, count, optional input SMILES file, uniqueness/randomization, and
Mol2Mol sampling strategy.

The `r4` conda environment exists and contains Python 3.10.17, REINVENT metadata 4.5.11, PyTorch
2.6.0, RDKit 2024.09.6, and pandas 2.2.3. Its `reinvent` console script is
`/Users/keetane/miniforge3/envs/r4/bin/reinvent`.

However, the environment's editable install maps `reinvent` and `reinvent_plugins` to the removed
path `/Users/keetane/Documents/apps/REINVENT4`. Running the current command therefore exits 1 with
`ModuleNotFoundError: No module named 'reinvent'`. Existing logs also record that older successful
runs used that old repository path.

With the live source explicitly supplied for diagnosis, an isolated two-molecule LibInvent smoke
run completed in 7.3 seconds with exit code 0. It wrote a valid two-row CSV and logged `Using
generator Libinvent`, `Sampling 2 SMILES`, and `Finished REINVENT`. This proves that the current
`libinvent.prior` and `child.smi` are mutually usable; it does not make the broken editable install
healthy.

## Inputs, outputs, and schemas

The live Streamlit runner uses these inputs:

- Mol2Mol: `input/parent.smi` and one of the Mol2Mol/Reinvent priors.
- LibInvent: `input/child.smi` and `priors/libinvent.prior`.
- LinkInvent: `input/warheads.smi` and `priors/linkinvent.prior`.

Observed raw output schemas differ by generator:

| Source | Columns |
| --- | --- |
| Mol2Mol raw | `SMILES`, `Input_SMILES`, `Tanimoto`, `NLL` |
| LibInvent raw smoke | `SMILES`, `Scaffold`, `R-groups`, `NLL` |
| Streamlit descriptor-enriched | raw columns plus `MW`, `LogP`, `HBD`, `HBA`, `TPSA`, `CtRotBonds`, `CtRings`, `CtAmides`, `CtAromaticRings`, `CtHetAromaticRings`, and sometimes `n_member_ring` |
| Analysis notebook output | `Rotatable_Bonds`, `Amide_Count`, `Ring_Count`, `Aromatic_Ring_Count`, and other analysis columns |

No inspected CSV contains `CtRotB`. The live writer creates `CtRotBonds`, while the viewer's sort
allowlist incorrectly names `CtRotB`; filtering names `CtRotBonds`. The future canonical field is
`RotB`, with aliases `RotB`, `RoB`, `CtRotB`, `CtRotBonds`, and the observed
`Rotatable_Bonds`. Conflicting simultaneous aliases must produce a warning rather than silently
choosing a value.

The legacy `Tanimoto` raw column and the Streamlit-added `Tanimoto_Similarity_to_Parent` are not
the same named contract. The latter is calculated with `Chem.RDKFingerprint` and
`DataStructs.TanimotoSimilarity`.

## Why `libinvent.csv` is absent

The latest `sampling.toml` targets `results/libinvent.csv` and was modified on 2026-10-16, but the
sampling log was last modified in 2025 and describes a successful Mol2Mol run from the removed old
repository path. A direct reproduction of the current subprocess command exits before logging or
sampling because its editable package target is missing. Then `run_reinvent()` ignores the nonzero
return code and immediately calls `pandas.read_csv(output_csv)`, converting the primary import
failure into a secondary file-not-found error. Other currently unguarded no-output paths include a
missing/invalid input fragment, missing prior, unsupported device, invalid TOML values, write
failure, or any nonzero REINVENT exit; all are masked by the same unconditional CSV read.

## Current risks and live-spec differences

- Sampling runs inside Streamlit and blocks the request/session process.
- Exit code is discarded; success is announced after only a CSV read and rewrite.
- One shared TOML and log path can be overwritten by concurrent sessions.
- Comment text is concatenated into filenames without traversal or collision protection.
- User-controlled subfolder input is joined directly to a server path.
- The full parent environment is inherited by the subprocess.
- Output is rewritten in place without atomic validation.
- Invalid SMILES and absent columns are incompletely handled.
- The repository's `environment.yml` lists only Python, numpy, and pandas; it does not reproduce
  the actual application or REINVENT environment.
- Documentation names `pages/auto_R_group.py`, but the live page is `pages/auto_R.py`.
- Existing outputs and TOMLs retain absolute paths under the former repository location.

## Phase 1 file plan

- `webapp/backend/app`: FastAPI entrypoint, configuration, error/API schemas, and service boundaries.
- `webapp/backend/app/db`: numbered SQLite migrations and per-operation connections with foreign
  keys, WAL, and busy timeout.
- `webapp/backend/app/models/job.py`: explicit status transition rules.
- `webapp/backend/worker`: independently launched worker process skeleton.
- `webapp/backend/tests`: migrations, health/OpenAPI, worker startup, and exhaustive transition tests.
- `webapp/frontend`: React/Vite shell with Vitest and Playwright foundations.

Phase 2 will add job persistence, atomic claim, cancellation, reconciliation, safe job directories,
and fake-command integration tests. Phase 3 will repair/package the REINVENT runtime contract and
implement the subprocess adapter using the command and output facts above.

## Decisions deferred to later phases

- Authentication/ownership policy is required before any non-loopback exposure.
- Tanimoto fingerprint parameters for newly computed similarity require explicit approval.
- Retention/deletion policy and deployment topology remain intentionally unset.
- Whether to repair `r4` in place or create a clean locked runtime is deferred to Phase 3; a clean
  runtime is recommended because the current editable installation is path-stale.

