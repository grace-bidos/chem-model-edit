# Worker Preflight Checker for Apptainer/MPI Single-Node

This document defines how to run the worker preflight checker introduced for GRA-72.

## Purpose

The checker verifies safe, non-destructive host prerequisites for v1 single-node worker execution:

- kernel and user namespace basics
- FUSE and overlay signals relevant to Apptainer
- Apptainer command availability and basic executable behavior
- environment indicators that imply unsupported multinode MPI assumptions

## Script

- Path: `scripts/preflight-worker-single-node.sh`
- Exit codes:
  - `0`: no failing checks (pass/warn only)
  - `1`: one or more failing checks

## Usage

Human-readable output:

```bash
scripts/preflight-worker-single-node.sh
```

Machine-readable JSON output:

```bash
scripts/preflight-worker-single-node.sh --json
```

Show help:

```bash
scripts/preflight-worker-single-node.sh --help
```

## Check Groups

### Kernel and namespace checks

- kernel release detection (`uname -r`)
- `/proc/sys/user/max_user_namespaces`
- `/proc/sys/kernel/unprivileged_userns_clone` (when exposed)

### FUSE and overlay checks

- `/dev/fuse` existence
- `fusermount3` or `fusermount` command presence
- overlay filesystem signal from `/proc/filesystems`

### Apptainer checks

- `apptainer` command availability
- `apptainer --version`
- `apptainer exec --help` probe

### MPI single-node scope checks

- Warn when multinode indicators are detected in environment variables such as:
  - `SLURM_NNODES`
  - `SLURM_JOB_NUM_NODES`
  - `OMPI_MCA_orte_num_nodes`
  - `HYDRA_NODELIST`
- This is expected for v1 behavior: multinode MPI is out of scope.

## JSON Output Contract

The JSON mode returns a single object with:

- `checker`: script identifier
- `version`: checker version string
- `timestamp_utc`: ISO-8601 timestamp
- `overall_status`: `pass`, `warn`, or `fail`
- `counts`: pass/warn/fail totals
- `checks`: list of per-check results:
  - `id`
  - `status`
  - `message`
  - `details`

## Safety Notes

- All probes are read-only or help/version command checks.
- The checker does not modify host settings, services, or filesystems.
