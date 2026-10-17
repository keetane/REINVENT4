# Phase 4 through Phase 6 implementation record

Date: 2026-10-16 (Asia/Tokyo)

## Phase 4

Parent input now supports direct SMILES plus a lazy, local Ketcher standalone editor. The backend is
the validation authority and returns canonical isomeric SMILES, MolBlock, formula, and atom count.
Job history, persistent selection, three-second polling, bounded logs, elapsed time, phase-only
progress, cancellation, and completed-result navigation are available. Closing the browser does not
affect job execution.

## Phase 5

Completed CSVs are normalized through one query service. Required descriptors use existing finite
values first and otherwise fixed RDKit calculations. Rotatable-bond aliases use the documented
precedence and conflicts produce a warning. Inclusive AND filters, allowlisted sorting, pagination,
summary statistics, histograms, deterministic 500-point scatter sampling, structure SVGs, and
filtered CSV/SDF exports share the same row set. Invalid SMILES are counted and skipped individually.

## Phase 6

The implementation received focused boundary tests, full backend regression, frontend unit/build
and browser tests, dependency audit, a real two-molecule LibInvent smoke, and a changed-code security
review. Operations, backup/recovery, staged migration, rollback, security limits, and direct
dependency licenses are documented. The legacy Streamlit files remain unchanged.

Rollback for these phases is to stop the three new processes and use Streamlit; no legacy data
migration or deletion is required.
