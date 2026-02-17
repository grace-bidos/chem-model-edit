# Submit -> Event -> Projection Lifecycle Smoke (GRA-105)

This runbook defines a lightweight, reproducible smoke path for local/dev validation of the runtime lifecycle contract:

- submit handling
- execution event normalization
- projection payload dispatch contract

## Goal

- Provide one command that validates the submit -> event -> projection lifecycle contract.
- Capture logs/artifacts in a stable location for quick diagnosis.
- Keep the smoke path independent from real Slurm infrastructure.

## Preconditions

Required commands:

- `uv`

Repository assumptions:

- Run from a checked-out repository that includes `apps/api` and tests.
- Python dependencies for `apps/api` are already synced.

Explicit non-requirements:

- No Slurm daemon, `sinfo`, or `srun` is required.
- No AiiDA profile/bootstrap is required.

## Entrypoint

- `scripts/submit-event-projection-smoke.sh`

## Commands

Default run:

```bash
scripts/submit-event-projection-smoke.sh
```

Custom artifact location:

```bash
scripts/submit-event-projection-smoke.sh \
  --artifact-dir investigations/artifacts/gra-105
```

Custom run label:

```bash
scripts/submit-event-projection-smoke.sh --run-label dev-loop
```

## Expected Output

PASS signals:

- Script exits with `0`.
- Output includes `Smoke result: PASS`.
- Output includes `Logs and artifacts: <path>`.
- Main pytest step prints `[ok] pytest-submit-event-projection`.

FAIL signals:

- Script exits with `1`.
- Output includes `Smoke result: FAIL`.
- Main pytest step prints `[fail] pytest-submit-event-projection`.

## Artifacts

Default artifact root:

- `investigations/artifacts/gra-105/`

Per-run directory format:

- `<UTC timestamp>-<run-label>` (example: `20260216T140500Z-local`)

Generated files:

- `manifest.txt`: run metadata (branch, commit, path)
- `pytest-submit-event-projection.log`: targeted lifecycle test output

## Troubleshooting

If script fails before pytest starts:

- Verify `uv` is installed and available in `PATH`.

If pytest step fails:

- Inspect `pytest-submit-event-projection.log` in the run artifact directory.
- Re-run the exact command from `[run]` output to iterate locally.

If imports/settings fail in `apps/api`:

- Sync dependencies in API project, then rerun:

```bash
cd apps/api
uv sync
```

If OpenAPI-related contract assertions fail unexpectedly:

- Confirm generated contract artifacts are current for your branch before rerunning smoke.
