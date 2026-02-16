# AiiDA-on-Slurm Docker Smoke Path (GRA-84)

This runbook defines a minimal smoke path for validating AiiDA + Docker runtime wiring and then checking Slurm scheduler reachability from the same host.

## Goal

- Validate the Docker-backed AiiDA runtime bootstrap (`core.psql_dos` profile).
- Validate that a `core.direct` computer test succeeds as a baseline.
- Attempt `core.slurm` computer smoke and explicitly report blockers.
- Provide an executable fallback path when Slurm control plane is unavailable.

## Preconditions

Required commands:

- `uv`
- `docker`
- `srun`
- `sinfo`

Required files:

- Environment file for bootstrap, default path: `apps/api/.env.aiida.dev`
- If missing, copy from `apps/api/.env.aiida.dev.example` and adjust values.

Operational assumptions:

- Docker daemon is running.
- Slurm client is installed on the host where smoke runs.
- For full Slurm success, the host can resolve and reach a Slurm controller (`slurmctld`) via local `slurm.conf` or DNS SRV configuration.

## Script

- `scripts/aiida-slurm-smoke-docker.sh`

Modes:

- `docker-slurm` (default): full path with Docker bootstrap + Slurm attempt
- `local-only`: fallback path with local SQLite profile and `core.direct` smoke only

## Commands

Full path (expected to return `0` or `2`):

```bash
scripts/aiida-slurm-smoke-docker.sh \
  --mode docker-slurm \
  --env-file apps/api/.env.aiida.dev
```

Fallback path (expected `0` if local scheduler wiring is OK):

```bash
scripts/aiida-slurm-smoke-docker.sh --mode local-only
```

Optional artifact location override:

```bash
scripts/aiida-slurm-smoke-docker.sh \
  --mode docker-slurm \
  --env-file apps/api/.env.aiida.dev \
  --artifact-dir investigations/artifacts/gra-84
```

## Exit Codes

- `0`: smoke passed for selected mode
- `1`: hard failure (bootstrap or direct scheduler smoke failed)
- `2`: Docker + direct smoke passed, Slurm stage blocked/failed

## What the Script Executes

In `docker-slurm` mode:

1. Runs `scripts/bootstrap-aiida-runtime.sh --apply` (container/network bootstrap only).
2. Creates/reuses `core.psql_dos` profile with current `verdi profile setup` options.
3. Ensures `core.direct` computer (`gra84-local-direct`) and runs `verdi computer test`.
4. Ensures `core.slurm` computer (`gra84-local-slurm`) and runs:
   - `sinfo`
   - `verdi computer test gra84-local-slurm`
5. If Slurm stage fails, emits blocker text and fallback command.

In `local-only` mode:

1. Creates/reuses a local `core.sqlite_dos` profile (`gra84-local-smoke`).
2. Ensures `core.direct` computer (`gra84-local-direct`) and runs `verdi computer test`.

## Blocker Interpretation

A typical blocker for local laptops/workstations without cluster config:

- `sinfo: fatal: Could not establish a configuration source`

This indicates the Slurm client exists but cannot find usable controller configuration (`slurm.conf`/DNS SRV).

## Artifacts

Logs are written under:

- `investigations/artifacts/gra-84/<timestamp>-<mode>/`

Key logs:

- `bootstrap-aiida-runtime.log`
- `aiida-computer-test-direct.log`
- `slurm-sinfo.log`
- `aiida-computer-test-slurm.log`

## Next Actions After Blocked Slurm Stage

1. Confirm Slurm controller config path and DNS SRV setup.
2. Verify host-level `sinfo` succeeds.
3. Re-run `--mode docker-slurm` and expect exit code `0`.
