# GRA-17: Tenant Isolation and Audit Logging Boundaries

- Status: Draft (Ask)
- Linear: `GRA-17`
- Scope: Boundary contract only (no runtime code in this item)
- Related ADRs:
  - `docs/adr/ADR-0001-system-of-record-boundaries.md`
  - `docs/adr/ADR-0003-projection-boundary-and-single-write-cutover.md`

## Objective

Define mandatory tenant isolation and audit logging boundaries for backend refresh
work so future Show PRs can implement behavior without introducing ownership
ambiguity or cross-tenant leakage.

## ADR Alignment

- `ADR-0001`: `Convex` remains product-facing SoR, `AiiDA/Postgres` remains
  execution-internal SoR, FastAPI is adapter/orchestration, and Redis is
  non-durable transitional infrastructure.
- `ADR-0003`: no long-lived dual-write ambiguity and no Redis durable ownership.
  Audit logging must be durable and queryable without redefining SoR boundaries.

## In Scope

- Tenant-context propagation and enforcement boundaries across adapter and
  product projection surfaces.
- Audit event contract and ownership boundaries for write/read paths.
- Acceptance criteria for the subsequent Show-layer implementation PR.

## Out of Scope

- Runtime implementation details (schema migrations, endpoint code, infra code).
- Introducing a new ADR (this spec aligns to accepted ADRs, it does not supersede
  them).
- Redefining job-state semantics from `ADR-0002`.

## Boundary Rules (Normative)

### BR-001: Mandatory Tenant Context

- All authenticated product-surface requests MUST resolve exactly one
  `tenant_id`.
- Requests with missing, malformed, or ambiguous `tenant_id` MUST be rejected
  before business mutation logic.
- No fallback "default tenant" behavior is allowed.

### BR-002: Ownership and SoR Mapping

| Concern                                        | Durable owner                 | Write path class                     | Read path class                         |
| ---------------------------------------------- | ----------------------------- | ------------------------------------ | --------------------------------------- |
| Product projections and collaboration metadata | `Convex`                      | FastAPI adapter -> Convex            | Product APIs/UI via Convex-backed reads |
| Scientific execution internals and provenance  | `AiiDA/Postgres`              | Execution pipeline -> AiiDA/Postgres | Detail endpoints via adapter            |
| Audit trail of tenant-scoped product actions   | Dedicated durable audit store | FastAPI adapter -> audit store       | Audit query APIs/reports                |
| Cache/lease/transient delivery                 | Redis (transitional only)     | Worker/adapter transient writes      | Operational runtime only                |

Boundary constraint: Redis MUST NOT become durable owner for tenant identity,
authorization decisions, or audit history.

### BR-003: Tenant-Scoped Product Data

- New or updated product-facing records MUST include `tenant_id`.
- Product reads and writes MUST apply tenant filters by default.
- Cross-tenant product queries are prohibited for tenant-scoped callers.
- Any intentionally shared dataset MUST be explicitly documented with rationale,
  owner, and read-scope policy.

### BR-004: Audit Event Minimum Contract

Each mutation-attempt audit event MUST include at least:

- `event_id` (stable unique id for idempotency)
- `request_id` (propagated correlation id)
- `tenant_id`
- `actor_id`
- `operation` (create/update/delete/execute/etc.)
- `resource_type`
- `resource_id` (if known at emission time)
- `outcome` (`attempted`, `succeeded`, or `failed`)
- `occurred_at` (UTC timestamp)

Audit events MAY include redacted metadata payloads, but MUST NOT contain
secrets in clear text.

### BR-005: Audit Durability and Ordering

- For mutation paths, audit emission MUST occur before success is acknowledged to
  the caller.
- If audit persistence cannot be guaranteed, the mutation path MUST NOT report
  success.
- Audit writes MUST be idempotent for retried requests.

### BR-006: Audit Query Access Boundaries

- Tenant-scoped actors can only query their own tenant audit records.
- Cross-tenant audit queries are platform-operator only and require explicit
  privileged role checks.
- Product-facing APIs MUST NOT expose unrestricted cross-tenant audit search.

### BR-007: Show PR Evidence Requirements

The Show-layer PR implementing these boundaries MUST include:

- SoR mapping table in PR body (`field/aggregate`, `owner`, `write path`,
  `read path`) per `ADR-0001`.
- Explicit statement of audit store owner and retention policy.
- Proof that Redis is not used as durable owner for audit or tenant semantics.
- Rollback notes describing behavior if audit store availability degrades.

## Acceptance Criteria

- AC-001: This Ask spec is merged with explicit references to
  `ADR-0001` and `ADR-0003`.
- AC-002: Boundary rules `BR-001` through `BR-007` are present and unambiguous.
- AC-003: A single ownership mapping exists for tenant projections, execution
  internals, and audit trail (no dual durable ownership for the same concern).
- AC-004: Tenant propagation and rejection behavior are explicitly defined for
  missing/invalid tenant context.
- AC-005: Audit event required fields are specified, including idempotency and
  outcome semantics.
- AC-006: Cross-tenant audit access is restricted to explicit privileged
  operators; tenant-scoped callers cannot query other tenants.
- AC-007: Show PR evidence requirements are defined and are sufficient to review
  SoR compliance.

## Risks and Follow-up Notes

- If a future Show PR needs to change durable owner assignments, that change
  requires a new ADR, not an implementation-only deviation.
- If shared/global datasets are introduced later, they need explicit exception
  documentation to avoid accidental cross-tenant leakage.
