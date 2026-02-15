# ADR-0004: Runtime Contract Boundary for BYO AiiDA/Slurm

- Status: Accepted
- Date: 2026-02-15
- Linear: `GRA-64`
- Owners: Backend platform
- Related: `ADR-0001`, `ADR-0003`, `GRA-67`, `GRA-68`, `GRA-69`

## Context

Backend refresh now assumes a BYO compute topology where users register a management node that hosts `AiiDA + Slurm`, and worker nodes execute containerized workloads.

Without a strict runtime contract boundary, parallel implementation lanes (Convex expansion, FastAPI orchestration hardening, Redis/RQ retirement) will conflict and regress ownership rules.

## Decision

### System responsibilities

- `Convex` is Product SoR.
  - Owns user-visible projection state, tenant-scoped query models, and product-level audit context.
  - Does not own scientific execution internals.
- `FastAPI` is orchestration adapter.
  - Owns request validation, command acceptance/rejection, idempotency handling, and projection update emission.
  - Is not a durable SoR for business state.
- `Management node stack (AiiDA/Postgres + Slurm)` is Execution SoR.
  - Owns scheduler/runtime execution details, queue execution behavior, provenance graph, and execution artifacts metadata.
  - Scientific execution internals and provenance durability remain owned by `AiiDA/Postgres` consistently with `ADR-0001`.

### Runtime contract surfaces

- `SubmitJobCommand` (`Convex`/client-facing call -> `FastAPI`)
- `ExecutionEvent` (`management node` -> `FastAPI`)
- `ProjectionUpdate` (`FastAPI` -> `Convex`)

Filename mapping in `docs/contracts/`:

- `SubmitJobCommand`: `command-submit-job.md`
- `ExecutionEvent`: `event-execution-lifecycle.md`
- `ProjectionUpdate`: `projection-update.md`

Each surface must define: required fields, idempotency key, retry rule, tenant boundary, and error semantics.

Operational note:

- Initial and updated contract definitions must live in `docs/contracts/`.
- One Markdown file per contract surface is required.

### Framework stance

- Keep `FastAPI` as the runtime adapter framework for this phase.
- Do not switch to Django/Django Ninja during the contract-fix phase.
- Re-evaluate framework choice only if we introduce Django-centric operational domains that materially exceed API orchestration needs.

### Migration stance for Redis/RQ

- Redis/RQ must not remain owners of business semantics.
- Redis may remain only as short-lived cache/lease/transport during transition.
- RQ paths are retirement targets, not expansion targets.
- Until retirement completion, Redis lease usage is transitional only and does not change ownership boundaries.

## Consequences

- Parallel lanes can proceed without ownership ambiguity:
  - Convex projection lane
  - FastAPI adapter lane
  - Management-node integration lane
  - Redis/RQ retirement lane
- Contract changes become explicit review points rather than incidental code drift.
- Pre-production phase prioritizes clear contracts over backward-compatibility complexity.

## Verification

Implementation PRs that touch runtime flow must include:

- Which contract surface changed (`SubmitJobCommand`, `ExecutionEvent`, or `ProjectionUpdate`)
- SoR mapping for changed fields
- Idempotency and retry behavior
- Tenant/auth boundary impact

Reject PRs that:

- reintroduce Redis/RQ business ownership
- collapse Product SoR and Execution SoR
- change contract semantics without updating contract docs in `docs/contracts/`
