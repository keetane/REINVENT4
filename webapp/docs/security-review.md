# Security review and known limits

Date: 2026-10-16 (Asia/Tokyo)

## Trust boundary

The supported boundary is one trusted local user on loopback. The browser is untrusted input to the
API, and REINVENT output is treated as untrusted structured data. The worker has the local user's
filesystem and model-runtime privileges. There is no identity or ownership layer, so the API must
not be exposed beyond the local machine.

## Implemented controls

- Model selection and run options are server allowlists. Subprocesses use an argument vector with
  `shell=False`, a minimal environment, an explicit source root, and a UUID job working directory.
- User comments do not become paths directly. Names are normalized and bounded, writes are confined
  beneath the job directory, and finalization uses atomic replacement.
- Parent SMILES are parsed and canonicalized by RDKit; empty, reaction, multi-fragment, invalid, and
  over-limit structures are rejected before TOML generation.
- Result paths are resolved beneath the configured repository, symlinks are rejected, CSV size is
  limited to 100 MiB, rows to 100,000, and the `SMILES` header is required.
- Result sort/filter/scatter fields are allowlisted. Numeric bounds must be finite. The service does
  not evaluate client-provided expressions.
- API errors omit tracebacks and local paths. Log responses are byte/line bounded and redact known
  configured paths plus common secret assignments.
- React's normal escaping is retained for comments and values. RDKit SVG is percent-encoded into an
  image data URL rather than inserted as HTML.
- CSV export neutralizes nonnumeric cells beginning with `=`, `+`, `-`, or `@`. SDF properties use
  fixed descriptor names.
- CORS defaults only to the two explicit loopback development origins. Runtime and frontend
  dependencies are locked.

## Known limits

- No authentication, authorization, per-job ownership, CSRF protection, rate limiting, audit log,
  retention policy, or encryption-at-rest is present.
- Results are reparsed and descriptors recomputed per request. Limits bound the input, but a trusted
  local user can still create CPU-heavy repeated requests.
- The result parser supports up to 100 MiB/100,000 rows; operating memory can exceed file size.
- Raw on-disk subprocess logs are not redacted. Only the API view is sanitized.
- Ketcher 3.18's lazy browser bundle is large (about 28.8 MB minified, 8.4 MB gzip). Direct SMILES
  entry remains usable while it loads or if initialization fails.
- Ketcher 3.18 declares Node 24.14.1 or newer. The implementation passed local tests on Node 22.17,
  but fresh installs should use the declared Node version.
- Tanimoto is not newly calculated because fingerprint type and parameters have not been approved.
- Automatic retention and deletion are deliberately absent.

## Review before broader deployment

Add authentication and object-level authorization for every list/get/cancel/result endpoint; bind
jobs to an owner; define CSRF and cookie/token handling; terminate TLS; add request and job quotas;
define audit, backup, and retention policy; review reverse-proxy limits; and repeat dependency,
license, penetration, and threat-model review. Do not rely on CORS as access control.
