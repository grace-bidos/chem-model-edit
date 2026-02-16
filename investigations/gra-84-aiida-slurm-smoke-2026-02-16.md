# GRA-84 AiiDA-on-Slurm Docker Smoke Investigation (2026-02-16)

## Scope

- Lane D smoke-path scaffold for Docker-backed AiiDA + Slurm validation.
- Keep outcome executable even when full Slurm path is blocked.

## Preflight

Executed in `/home/grace/projects/chem-model-edit/.worktrees/gra-84`:

- `for tool in node corepack pnpm uv; do command -v "$tool" >/dev/null || exit 1; done` -> pass
- `node -v` -> `v22.16.0`
- `pnpm -v` -> `10.27.0`
- `uv --version` -> `uv 0.6.14`

Runtime availability checks:

- `command -v docker` -> `/usr/bin/docker`
- `command -v srun` -> `/usr/bin/srun`
- `command -v sbatch` -> `/usr/bin/sbatch`

## What Ran

### 1) Fallback smoke (local-only)

Command:

```bash
scripts/aiida-slurm-smoke-docker.sh --mode local-only
```

Result:

- Exit code: `0`
- `core.sqlite_dos` local profile path is usable.
- `verdi computer test` for `gra84-local-direct` passed (`6/6` tests).

Artifact examples:

- `investigations/artifacts/gra-84/20260216-090739-local-only/aiida-computer-test-direct.log`

### 2) Docker + Slurm smoke

Command:

```bash
scripts/aiida-slurm-smoke-docker.sh --mode docker-slurm --env-file apps/api/.env.aiida.dev.example
```

Result:

- Exit code: `2` (blocked at Slurm stage by design)
- Docker runtime bootstrap succeeded (PostgreSQL + RabbitMQ containers available).
- AiiDA `core.psql_dos` profile (`chem-model-dev`) was set up successfully.
- Direct scheduler test (`gra84-local-direct`) passed.
- Slurm stage failed at controller/config resolution.

Artifact examples:

- `investigations/artifacts/gra-84/20260216-090718-docker-slurm/bootstrap-aiida-runtime.log`
- `investigations/artifacts/gra-84/20260216-090718-docker-slurm/aiida-computer-test-direct.log`
- `investigations/artifacts/gra-84/20260216-090718-docker-slurm/slurm-sinfo.log`
- `investigations/artifacts/gra-84/20260216-090718-docker-slurm/aiida-computer-test-slurm.log`

## Blocker

Exact blocker observed:

- `sinfo: fatal: Could not establish a configuration source`
- `squeue: fatal: Could not establish a configuration source`

Interpretation:

- Slurm client binaries exist, but host-side Slurm control-plane configuration is missing or unreachable (`slurm.conf` / DNS SRV / controller reachability).

## Executable Fallback Path

Use local scheduler smoke first:

```bash
scripts/aiida-slurm-smoke-docker.sh --mode local-only
```

Then, after fixing host Slurm config/controller access, retry full path:

```bash
scripts/aiida-slurm-smoke-docker.sh --mode docker-slurm --env-file apps/api/.env.aiida.dev
```

## Notes

- Existing helper `scripts/bootstrap-aiida-runtime.sh --init-profile` currently uses outdated AiiDA CLI option names and fails profile initialization in this environment.
- The new smoke script works around this by:
  - using bootstrap only for Docker runtime setup,
  - then creating `core.psql_dos` profile with current option names.
