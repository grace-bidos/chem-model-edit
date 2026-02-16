# Managed AiiDA Runtime Health/Readiness Checks

This feature adds executable managed AiiDA runtime probes to FastAPI health endpoints.

## Endpoints

- `GET /api/health`
  - Always returns `200`.
  - `status` is:
    - `ok` when all enabled checks are `ok` (or `skipped`)
    - `degraded` when any enabled check fails
- `GET /api/ready`
  - Returns `200` when ready.
  - Returns `503` when any enabled readiness check fails.
  - `status` is `ready` or `not_ready`.

Both endpoints include:

- `checks.managed_aiida_runtime.status`: `ok | failed | skipped`
- `checks.managed_aiida_runtime.detail`
- `checks.user_managed_deep_readiness.status`: `ok | failed | skipped`
- `checks.user_managed_deep_readiness.detail`
- `checks.user_managed_deep_readiness.probes` with per-command probe results:
  - `scontrol_ping` (`scontrol ping`)
  - `sinfo` (`sinfo`)
  - `verdi_status` (`verdi status`)
  - `verdi_profile` (`verdi profile list`)
- Optional failure context:
  - `error` (`misconfigured`, `unreachable`, `upstream_unhealthy`)
  - `http_status`
  - `url`
  - `latency_ms`

Deep readiness probe failures additionally report:

- `error` (`tooling_missing`, `command_failed`, `timeout`, `misconfigured`)
- `failed_probe` (the first failing probe key)
- per-probe `command`, optional `exit_code`, `stdout`, `stderr`

## Environment Variables

- `AIIA_MANAGED_CHECKS_ENABLED`
  - `true/false` (default: disabled)
- `AIIA_MANAGED_BASE_URL`
  - Required when checks are enabled.
- `AIIA_MANAGED_HEALTH_PATH`
  - Default: `/health`
- `AIIA_MANAGED_READY_PATH`
  - Default: `/ready`
- `AIIA_MANAGED_TIMEOUT_SECONDS`
  - Default: `3`
- `AIIA_MANAGED_BEARER_TOKEN`
  - Optional bearer token for probe requests.
- `AIIA_USER_MANAGED_DEEP_READY_ENABLED`
  - `true/false` (default: disabled)
  - Opt-in deep readiness checks for user-managed AiiDA+Slurm runtime.
- `AIIA_USER_MANAGED_DEEP_READY_TIMEOUT_SECONDS`
  - Default: `5`
  - Per-command timeout for deep readiness probes.
- `AIIDA_PROFILE`
  - Optional profile name to verify in `verdi profile list`.
  - When set and missing from output, deep readiness fails with `misconfigured`.

## Quick Validation

Disabled mode (default):

```bash
curl -s http://localhost:8000/api/health
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8000/api/ready
```

Expected:

- `/api/health` => `status: ok`, check status `skipped`
- `/api/ready` => HTTP `200`

Deep readiness enabled and missing Slurm tooling:

```bash
export AIIA_USER_MANAGED_DEEP_READY_ENABLED=true
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8000/api/ready
```

Expected:

- If required commands are missing (for example `scontrol`), `/api/ready` => HTTP `503`
- `checks.user_managed_deep_readiness.error` => `tooling_missing`

Enabled + missing base URL:

```bash
export AIIA_MANAGED_CHECKS_ENABLED=true
unset AIIA_MANAGED_BASE_URL
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8000/api/ready
```

Expected:

- `/api/ready` => HTTP `503`
- check status `failed` with `error: misconfigured`

Enabled + healthy upstream:

```bash
export AIIA_MANAGED_CHECKS_ENABLED=true
export AIIA_MANAGED_BASE_URL=http://<managed-aiida-host>
curl -s http://localhost:8000/api/ready
```

Expected:

- HTTP `200`
- check status `ok`
