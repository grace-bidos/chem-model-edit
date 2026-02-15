# ADR-0001: System-of-Record Boundaries and Ownership Model

- Status: Accepted
- Date: 2026-02-15
- Linear: `GRA-12`
- Owners: Backend platform

## Context

Backend refresh introduces multiple persistence and orchestration surfaces (`Convex`, `AiiDA/Postgres`, FastAPI adapter, and transitional Redis paths). Without explicit ownership, runtime behavior can diverge across systems and create ambiguous truth sources.

## Decision

Define clear ownership boundaries:

- `AiiDA/Postgres` is the system of record for scientific workflow execution internals.
- `Convex` is the system of record for product-facing job projection and UI read models.
- `FastAPI` is a thin adapter/orchestration layer and is not a long-term source of truth.
- `Redis` is transitional transport/cache/lease infrastructure only; it must not be treated as durable truth.

Ownership by concern:

- Job execution internals (scheduler/runtime metadata): `AiiDA/Postgres`
- User-visible job lifecycle and query model: `Convex`
- Authentication/session data (current implementation): dedicated auth store
- Queue lease and transient delivery state: Redis (temporary)

## Consequences

- All new backend-refresh features must choose one durable owner (`AiiDA/Postgres` or `Convex`) for each persisted field.
- Adapter-level writes must be projection or command-routing oriented, not ownership-oriented.
- Migration work (`GRA-21`, `GRA-22`) must remove Redis-owned business semantics.

## Implementation Notes

- `GRA-14`, `GRA-16`, `GRA-18`, `GRA-20` must preserve the ownership split above.
- API contracts exposed to clients should be modeled from Convex-facing read models.
- Scientific execution write paths should terminate in AiiDA/Postgres-owned structures.

## Verification

- For each runtime PR, include a short "SoR mapping" table in PR description:
  - field/aggregate
  - owner
  - write path
  - read path
- Reject PRs that introduce dual-write ambiguity without an explicit migration contract.
