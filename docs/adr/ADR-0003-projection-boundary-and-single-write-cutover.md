# ADR-0003: Projection Boundary and Single-Write Cutover Policy

- Status: Accepted
- Date: 2026-02-15
- Linear: `GRA-37`
- Owners: Backend platform
- Related: `GRA-16`, `GRA-20`, `GRA-21`, `GRA-22`, `GRA-24`

## Context

Backend refresh must keep product UX fast while preserving AiiDA strengths (provenance-rich execution internals). If we mirror AiiDA internals broadly into Convex, migration and feature cost will grow with every workflow detail added. We also need a clear migration policy from Redis that avoids prolonged dual-write ambiguity.

## Decision

Projection and ownership policy:

- `Convex` remains product-facing system of record for UI read models and collaboration-facing metadata.
- `AiiDA/Postgres` remains system of record for scientific execution internals and provenance graph details.
- Use stable cross-reference keys (for example `aiida_node_id`) to join product projections with execution details.
- Deep execution detail views must read from AiiDA through adapter endpoints on demand.
- Do not bulk-mirror full AiiDA execution/provenance models into Convex.

Field-level ownership clarifications for pending work:

- Queue target profiles and active profile selection: `Convex` owner.
- Worker authentication credentials:
  - durable owner: dedicated auth store
  - token format: opaque tokens
  - storage: hashed-at-rest representation
  - Redis may be used only as short-lived cache, not durable authority.

Cutover policy:

- Use single-write migration policy for runtime paths (`GRA-21`, `GRA-22`); no long-lived dual-write phase.
- Transitional emergency-stop controls are allowed only as short-lived operations switches with explicit owner and expiry in PR description.
- Redis stays limited to lease/transport/transient cache semantics.

## Consequences

- Product list/status flows stay fast by reading Convex projection models.
- Detailed execution/provenance views require adapter calls to AiiDA sources.
- Model-drift risk is reduced because Convex projection scope is intentionally narrow.
- Migration sequencing becomes stricter: ownership must be explicit before write-path cutover.

## Implementation Notes

- `GRA-16`: event relay writes only product projection fields to Convex.
- `GRA-20`: adapter surfaces summary-vs-detail boundaries explicitly.
- `GRA-21`: cutover flags include temporary emergency-stop fields with expiry metadata.
- `GRA-22`: remove Redis-owned business semantics from result/meta/ownership/auth paths.
- `GRA-24`: restart durability must come from durable owners (`Convex` and/or `AiiDA/Postgres`), not Redis keys.

## Verification

- Runtime PRs must include:
  - projection field list (what is stored in Convex)
  - detail read path (where AiiDA is queried)
  - cutover plan (single-write switch step and rollback behavior)
- Reject PRs that:
  - introduce broad AiiDA-to-Convex model mirroring
  - reintroduce Redis as durable owner for business semantics
  - add emergency-stop flags without explicit expiry/removal plan.
