# Redis Worker Queue/Result Retirement Plan (`GRA-22`)

This document defines the concrete retirement scope for Redis-owned ZPE worker
queue/result paths and the final removal checklist.

It follows ADR-0001 and ADR-0003:

- no long-lived dual-write path
- no Redis ownership of durable business semantics
- keep rollback bounded to short-lived operational controls

## Scope

`GRA-22` covers retirement planning and removal sequencing for legacy Redis
paths used by:

- job enqueue/status/result read
- worker lease/result/failure callbacks
- worker enrollment/auth and queue target selection
- temporary cutover flags currently stored in Redis

Actual endpoint/service deletion should happen only after next-gen
submission/lease/result/ownership paths are production-ready.

## Deletion and Deprecation Inventory

### Router-level legacy endpoints

- `apps/api/app/routers/zpe.py`
  - Legacy user job endpoints:
    - `POST /api/zpe/jobs`
    - `GET /api/zpe/jobs/{job_id}`
    - `GET /api/zpe/jobs/{job_id}/result`
    - `GET /api/zpe/jobs/{job_id}/files`
  - Legacy worker endpoints:
    - `POST /api/zpe/compute/enroll-tokens`
    - `POST /api/zpe/compute/servers`
    - `POST /api/zpe/compute/jobs/lease`
    - `POST /api/zpe/compute/jobs/{job_id}/result`
    - `POST /api/zpe/compute/jobs/{job_id}/failed`
  - Legacy queue target endpoints:
    - `GET /api/zpe/targets`
    - `PUT /api/zpe/targets/{target_id}/active`
  - Legacy migration control endpoints:
    - `GET /api/zpe/admin/ops`
    - `PATCH /api/zpe/admin/ops`
  - Keep as temporary exception:
    - `DELETE /api/zpe/compute/servers/{server_id}` for emergency revoke during
      transition; remove once replacement auth/revoke flow is active.

### Redis-backed service modules

- `apps/api/services/zpe/queue.py`
- `apps/api/services/zpe/backends.py`
- `apps/api/services/zpe/http_queue.py`
- `apps/api/services/zpe/lease.py`
- `apps/api/services/zpe/compute_results.py`
- `apps/api/services/zpe/result_store.py`
- `apps/api/services/zpe/job_owner.py`
- `apps/api/services/zpe/job_meta.py`
- `apps/api/services/zpe/enroll.py`
- `apps/api/services/zpe/worker_auth.py`
- `apps/api/services/zpe/queue_targets.py`
- `apps/api/services/zpe/ops_flags.py`

### Settings and exports to retire

- `apps/api/services/zpe/settings.py`
  - retire Redis-only knobs once replacements land:
    - `redis_url`
    - `result_store=redis`
    - legacy cutover defaults:
      - `cutover_submission_route`
      - `cutover_result_read_source`
      - `legacy_worker_endpoints_enabled`
- `apps/api/services/zpe/__init__.py`
  - remove exports for deleted legacy modules

### Docs/tests to update with final removal

- `docs/process/redis-worker-cutover-flags.md`
- `docs/zpe-worker-setup.md`
- Redis-focused tests under `apps/api/tests/test_zpe_*`
- generated API contract (`packages/api-client/src/generated/schema.ts`)

## Cleanup Migration Notes

### Required order

1. Replace write/read paths first.
   - next-gen submission path replaces `POST /api/zpe/jobs`
   - projection/AiiDA-backed read path replaces `GET /api/zpe/jobs/*`
2. Replace worker runtime contract.
   - next-gen lease/result/failure flow replaces legacy `/compute/jobs/*`
3. Replace auth and ownership stores.
   - move worker auth/enrollment and job owner/meta to durable authorities
4. Disable legacy routes through cutover controls.
   - set:
     - `submission_route=next-gen`
     - `result_read_source=projection`
     - `legacy_worker_endpoints_enabled=false`
5. Remove Redis-bound modules and endpoint wiring in one cleanup PR layer.

### Rollback policy

- Rollback is operational only and short-lived:
  - temporarily re-enable legacy route controls if replacement path regresses
- do not add new dual-write behavior
- keep emergency revoke available until replacement credential lifecycle proves
  stable

### Out-of-scope for this planning PR

- introducing replacement next-gen services
- changing production defaults before replacement path verification
- deleting legacy code before readiness gates are met

## Final Removal PR Checklist

- [ ] Router cutover complete: `apps/api/app/routers/zpe.py` no longer imports or
      references legacy Redis worker modules.
- [ ] Legacy endpoints removed (or replaced) for `/api/zpe/jobs*`,
      `/api/zpe/compute/jobs*`, `/api/zpe/targets*`, and `/api/zpe/admin/ops`.
- [ ] Redis-backed service modules deleted:
      `queue.py`, `backends.py`, `http_queue.py`, `lease.py`,
      `compute_results.py`, `result_store.py`, `job_owner.py`, `job_meta.py`,
      `enroll.py`, `worker_auth.py`, `queue_targets.py`, `ops_flags.py`.
- [ ] `apps/api/services/zpe/__init__.py` exports updated for post-retirement API.
- [ ] `apps/api/services/zpe/settings.py` Redis-only and legacy-cutover settings
      removed or replaced.
- [ ] API contract regenerated and committed (`packages/api-client/openapi/openapi.json`,
      `packages/api-client/src/generated/schema.ts`).
- [ ] Legacy setup/operations docs updated to next-gen runtime instructions.
- [ ] Redis dependency audit passes for retired domains:
      `rg -n \"redis|Redis\" apps/api/services/zpe apps/api/app/routers/zpe.py`.
- [ ] Validation passes:
      `uv run --project apps/api pytest`,
      `uv run --project apps/api ruff check apps/api`,
      `uv run --project apps/api mypy apps/api`.
