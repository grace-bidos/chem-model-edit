# GRA-72 Worker Preflight Checker Notes

Date: 2026-02-16
Branch: `feature/gra-72-worker-preflight-single-node`

## Objective

Ship a safe, host-side preflight checker for v1 worker prerequisites with dual output modes (human + JSON).

## Implementation Summary

Added:

- `scripts/preflight-worker-single-node.sh`
- `docs/process/worker-preflight-single-node.md`

Checker behavior:

- validates kernel/user namespace/FUSE/overlay signals
- validates Apptainer availability and `exec` subcommand response
- warns on multinode MPI assumptions in current environment
- supports `--json` for machine-readable automation

## Verification

Syntax check:

```bash
bash -n scripts/preflight-worker-single-node.sh
```

Runtime check (human):

```bash
scripts/preflight-worker-single-node.sh
```

Runtime check (json):

```bash
scripts/preflight-worker-single-node.sh --json
```

## Observed Output Summary (this lane)

- Overall status: `fail`
- Pass checks included:
  - kernel release detection
  - user namespaces enabled (`max_user_namespaces`)
  - `/dev/fuse` present
  - `fusermount3` command present
  - overlay filesystem detected
  - multinode MPI assumption not detected
  - `mpirun` command present
- Warning checks included:
  - `unprivileged_userns_clone` path not exposed
  - Apptainer exec probe skipped because Apptainer is missing
- Failing checks included:
  - `apptainer` command missing

Interpretation:

- Host appears close to ready but cannot run containerized worker tasks until Apptainer is installed.
- Single-node scope warning behavior is active and non-blocking.
