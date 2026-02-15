# Managed AiiDA Runtime Health/Readiness Checks

This feature adds executable managed AiiDA runtime probes to FastAPI health endpoints.

## Endpoints

- `GET /api/health`
  - Always returns `200`.
  - `status` is:
    - `ok` when managed AiiDA check is `ok` or `skipped`
    - `degraded` when managed AiiDA check fails
- `GET /api/ready`
  - Returns `200` when ready.
  - Returns `503` when managed AiiDA readiness probe fails.
  - `status` is `ready` or `not_ready`.

Both endpoints include:

- `checks.managed_aiida_runtime.status`: `ok | failed | skipped`
- `checks.managed_aiida_runtime.detail`
- Optional failure context:
  - `error` (`misconfigured`, `unreachable`, `upstream_unhealthy`)
  - `http_status`
  - `url`
  - `latency_ms`

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

## Quick Validation

Disabled mode (default):

```bash
curl -s http://localhost:8000/api/health
curl -s -o /dev/null -w "%{http_code}\n" http://localhost:8000/api/ready
```

Expected:

- `/api/health` => `status: ok`, check status `skipped`
- `/api/ready` => HTTP `200`

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
