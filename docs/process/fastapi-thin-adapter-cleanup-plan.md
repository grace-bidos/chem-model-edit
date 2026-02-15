# FastAPI Thin Adapter Cleanup Plan (`GRA-20`)

## Goal

Rescope FastAPI into a thin adapter layer for backend refresh while preserving domain-heavy logic in Python service modules and external systems of record.

This plan is aligned with:

- `docs/adr/ADR-0001-system-of-record-boundaries.md`
- `docs/adr/ADR-0003-projection-boundary-and-single-write-cutover.md`

## Deliverables Covered

- Route/service cleanup plan
- Retained responsibilities for Python domain wrappers
- Deprecated endpoint list (planned, not immediate removal)

## Target Layer Boundaries

FastAPI adapter must keep only:

- authentication/authorization gates
- request/response schema validation
- orchestration between application services
- error/status-code mapping
- request-scoped logging metadata wiring

FastAPI adapter must not own:

- durable business state semantics
- queue or lease internals
- worker token lifecycle internals
- scientific workflow provenance internals

## Retained Python Responsibilities

The following responsibilities remain in Python service modules and are intentionally not absorbed into route handlers:

- QE input parsing and normalization wrappers:
  - `apps/api/services/zpe/parse.py`
  - `apps/api/services/parse.py`
- ASE/pymatgen conversion and structure operations:
  - `apps/api/services/structures.py`
  - `apps/api/services/supercells.py`
  - `apps/api/services/transplant.py`
- QE runtime wrappers and execution utilities:
  - `apps/api/services/zpe/qe.py`
  - `apps/api/services/zpe/worker.py`
  - `apps/api/services/zpe/thermo.py`
- Result/lease/queue orchestration internals (transitional until cutover completes):
  - `apps/api/services/zpe/result_store.py`
  - `apps/api/services/zpe/lease.py`
  - `apps/api/services/zpe/compute_results.py`

## Route/Service Cleanup Plan

### Phase 1: Router Slimming (No Behavior Change)

- Extract per-endpoint application-level orchestration from `apps/api/app/routers/zpe.py` into dedicated application services under `apps/api/services/zpe/`.
- Keep route handlers focused on auth gates, DTO conversion, and HTTP error mapping.
- Keep current API contract stable.

### Phase 2: Ownership-Safe Refactor

- Ensure each persisted field has a single durable owner (`Convex` or `AiiDA/Postgres`).
- Prevent reintroduction of Redis durable ownership semantics.
- Keep Redis limited to transient lease/transport/cache paths.

### Phase 3: Endpoint Deprecation and Removal

- Mark selected legacy endpoints as deprecated in OpenAPI only after replacement paths are available.
- Remove deprecated endpoints after migration window and compatibility checks.

## Planned Deprecated Endpoints

These are planned deprecations, not immediate removals.

1. `GET /api/zpe/targets`
   - Reason: queue target profile reads move to product-facing projection owner.
   - Replacement direction: Convex-backed projection query.
   - Earliest removal window: after target read path migration in backend refresh runtime stack.

2. `PUT /api/zpe/targets/{target_id}/active`
   - Reason: active target selection ownership moves to product-facing projection owner.
   - Replacement direction: Convex-backed projection mutation.
   - Earliest removal window: after frontend/BFF migration to replacement mutation.

3. `POST /api/zpe/compute/enroll-tokens`
4. `POST /api/zpe/compute/servers`
5. `DELETE /api/zpe/compute/servers/{server_id}`
   - Reason: worker credential lifecycle should be owned by dedicated auth authority, not route-local Redis semantics.
   - Replacement direction: dedicated auth/worker enrollment service contract.
   - Earliest removal window: after token issuance/revocation migration is complete.

6. `GET /api/zpe/admin/ops`
7. `PATCH /api/zpe/admin/ops`
   - Reason: transitional emergency controls should move to managed platform flag ownership.
   - Replacement direction: platform-managed operational control plane.
   - Earliest removal window: after single-write cutover stabilization period.

## Validation Checklist for `GRA-20` PR Layers

- Route handlers do not contain durable state ownership logic.
- Route handlers do not implement queue/lease/token internals directly.
- Service extraction preserves current API responses and status codes.
- Deprecation targets are documented with migration destination and removal condition.
