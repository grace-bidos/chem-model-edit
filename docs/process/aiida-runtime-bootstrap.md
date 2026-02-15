# AiiDA Runtime Bootstrap (Dev/Staging)

This document defines a minimal bootstrap path for local dev and staging-like validation of the AiiDA runtime baseline used in backend refresh work.

Scope of this slice (`GRA-14`):

- Add optional AiiDA dependency group in `apps/api/pyproject.toml`.
- Provide a dedicated env template: `apps/api/.env.aiida.dev.example`.
- Provide a safe bootstrap helper: `scripts/bootstrap-aiida-runtime.sh`.

## Goals

- Keep setup repeatable across contributors.
- Default to non-destructive behavior.
- Make runtime prerequisites explicit (PostgreSQL + RabbitMQ + AiiDA profile config path).

## Prerequisites

- `uv` (Python dependency management)
- `docker` (for local PostgreSQL and RabbitMQ containers)

## Environment File

1. Copy the template to a local env file:
   - `cp apps/api/.env.aiida.dev.example apps/api/.env.aiida.dev`
2. Adjust values for your environment:
   - host/port/credentials
   - profile name
   - local runtime paths

## Bootstrap Script

The script is dry-run by default:

```bash
scripts/bootstrap-aiida-runtime.sh --env-file apps/api/.env.aiida.dev
```

Apply mode executes commands:

```bash
scripts/bootstrap-aiida-runtime.sh --apply --env-file apps/api/.env.aiida.dev
```

One-time profile setup can be added with:

```bash
scripts/bootstrap-aiida-runtime.sh --apply --env-file apps/api/.env.aiida.dev --init-profile
```

## Safety Model

- No destructive operations are performed.
- Existing Docker network/containers are reused.
- Existing containers are started (not recreated) when stopped.
- New resources are created only when missing.
- `--apply` is required for any command execution.

## Verification

After bootstrap:

```bash
uv sync --project apps/api --group aiida
AIIDA_PATH=.just-runtime/aiida/config uv run --project apps/api --group aiida verdi status
```

If `--init-profile` is used, the script attempts non-interactive `verdi profile setup core.psql_dos` only when the target profile is absent.
