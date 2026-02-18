# Modal + AiiDA/Slurm Runtime Gate

This runbook defines executable checks for the runtime-only compute contract (`/api/runtime/*`) on FastAPI deployed to Modal.

## Scope

- Runtime submit/status/detail/projection endpoints on Modal
- Tenant-scoped auth and proxy-auth behavior
- User-managed AiiDA/Slurm execution event flow through stateless runtime gateway
- Convex projection consistency checks

## Local Operator Setup

For WSL + Hyper-V local operator stability (secret file placement, env sync flow,
and VM networking pitfalls), see:

- `docs/process/local-dev-runtime-ops.md`
- `scripts/dev/bootstrap_dev_secrets_file.sh`
- `scripts/dev/apply_dev_runtime_env.sh`
- `scripts/dev/runtime_submit_smoke.sh`

## Preconditions

- Modal app deployed from `apps/api/modal_app.py`
- Proxy auth is enabled (`@modal.asgi_app(requires_proxy_auth=True)`)
- Cloudflare Worker and web BFF are configured with:
  - `API_BASE` (Modal endpoint)
  - `MODAL_PROXY_SECRET` (for proxy auth)
- Management node can post runtime events with valid auth + tenant headers

## Gate 1: Modal API reachability

```bash
curl -fsS https://<workspace>--chem-model-edit-api-api.modal.run/api/health
```

Expected:

- HTTP `200`
- `checks.managed_aiida_runtime` appears in payload when deep checks are enabled

## Gate 2: Runtime submit contract

```bash
curl -fsS -X POST \
  -H "Authorization: Bearer <JWT>" \
  -H "x-tenant-id: <tenant_id>" \
  -H "content-type: application/json" \
  https://<workspace>--chem-model-edit-api-api.modal.run/api/runtime/jobs:submit \
  -d '{
    "tenant_id": "<tenant_id>",
    "workspace_id": "<tenant_id>",
    "job_id": "job-gate-001",
    "idempotency_key": "submit-job-gate-001",
    "management_node_id": "mgmt-node-1",
    "execution_profile": {"queue_name": "default"},
    "resource_shape": {"cpu": 1, "memory_mib": 2048, "walltime_seconds": 3600},
    "payload_ref": {"input_uri": "inline://gate/job-gate-001"},
    "requested_by": {"user_id": "<jwt-sub>"}
  }'
```

Expected:

- HTTP `202`
- response includes `job_id`, `submission_id`, `accepted_at`, `trace_id`

## Gate 3: Runtime event ingestion

Post a running and completed event from management-node integration.

```bash
curl -fsS -X POST \
  -H "Authorization: Bearer <JWT>" \
  -H "x-tenant-id: <tenant_id>" \
  -H "content-type: application/json" \
  https://<workspace>--chem-model-edit-api-api.modal.run/api/runtime/jobs/job-gate-001/events \
  -d '{
    "event_id": "evt-running-001",
    "tenant_id": "<tenant_id>",
    "workspace_id": "<tenant_id>",
    "job_id": "job-gate-001",
    "submission_id": "<submission_id>",
    "execution_id": "aiida-<node-id>",
    "state": "running",
    "occurred_at": "2026-02-18T00:00:00Z",
    "trace_id": "trace-running-001"
  }'
```

Expected:

- HTTP `200`
- payload `{"ok": true, "idempotent": false}`

Repeat with `state: completed` and `result_ref.output_uri`.

## Gate 4: Runtime read checks

```bash
curl -fsS \
  -H "Authorization: Bearer <JWT>" \
  -H "x-tenant-id: <tenant_id>" \
  https://<workspace>--chem-model-edit-api-api.modal.run/api/runtime/jobs/job-gate-001

curl -fsS \
  -H "Authorization: Bearer <JWT>" \
  -H "x-tenant-id: <tenant_id>" \
  https://<workspace>--chem-model-edit-api-api.modal.run/api/runtime/jobs/job-gate-001/detail

curl -fsS \
  -H "Authorization: Bearer <JWT>" \
  -H "x-tenant-id: <tenant_id>" \
  https://<workspace>--chem-model-edit-api-api.modal.run/api/runtime/jobs/job-gate-001/projection
```

Expected:

- `/api/runtime/jobs/{id}` returns terminal state (`completed` or `failed`)
- `/detail` contains event history in order
- `/projection` returns app-facing status (`finished` when completed)

## Gate 5: Tenant isolation

- Re-run Gate 4 with a different tenant header/JWT.
- Expected: `404` or scope error, never cross-tenant data.

## Gate 6: Convex projection alignment

- Validate projected job status in Convex read model for `job-gate-001`.
- Expected: Convex status monotonicity matches runtime event sequence.

## Failure diagnostics

- `401/403` on runtime endpoints: JWT claims mismatch (`tenant_id`, `sub`) or missing proxy auth.
- `409` on events: idempotency conflict, duplicated `event_id` with payload drift.
- `404` after submit: wrong tenant scope during read.
- projection mismatch: inspect runtime event ordering and Convex relay logs.

## Retirement note

Legacy `/api/zpe/*` endpoints are intentionally removed. Runtime-facing entrypoints live under `/api/runtime/*`, and onboarding dry-run lives under `/api/onboarding/dry-run`. Any operational script relying on RQ workers is out of scope for runtime operation as of this cutover.
