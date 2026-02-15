# Redis Worker Cutover Flags (`GRA-21`)

This document defines staged migration controls for retiring Redis worker paths safely without long-lived dual-write behavior.

## Flags

Current operational flags exposed at `GET/PATCH /api/zpe/admin/ops`:

- `submission_enabled: bool`
- `dequeue_enabled: bool`
- `submission_route: "redis-worker" | "next-gen"`
- `result_read_source: "redis" | "projection"`
- `legacy_worker_endpoints_enabled: bool`

## Defaults

- `submission_enabled=true`
- `dequeue_enabled=true`
- `submission_route="redis-worker"`
- `result_read_source="redis"`
- `legacy_worker_endpoints_enabled=true`

Defaults preserve current behavior and are backward compatible.

## Staged Migration Policy

1. **Preparation**
   - Keep defaults.
   - Implement next-gen write/read paths in parallel (without enabling them).

2. **Submission Cutover**
   - Set `submission_route="next-gen"` only when next-gen submission path is production-ready.
   - Keep `result_read_source="redis"` until downstream read path is validated.

3. **Read Cutover**
   - Set `result_read_source="projection"` when projection reads are validated.

4. **Legacy Endpoint Disablement**
   - Set `legacy_worker_endpoints_enabled=false` after next-gen worker/ingest paths are active.
   - This disables legacy compute enrollment/lease/result/failure endpoints.
   - Admin revoke (`DELETE /api/zpe/compute/servers/{server_id}`) remains available for emergency token invalidation.

5. **Cleanup**
   - Remove disabled legacy endpoints in `GRA-22` once migration verification is complete.

## Guardrails

- `submission_route != "redis-worker"` currently blocks legacy `/api/zpe/jobs` enqueue path.
- `result_read_source != "redis"` currently blocks legacy `/api/zpe/jobs/*` status/result/file read paths.
- `legacy_worker_endpoints_enabled=false` blocks legacy compute enrollment/lease/result/failure endpoints.
- Admin revoke endpoint remains enabled to preserve incident response ability during and after cutover.
- Flags are intended as short-lived migration controls, not permanent configuration.
