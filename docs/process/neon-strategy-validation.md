# Neon Strategy Validation for Durable Data (`GRA-19`)

This document validates whether Neon should be adopted as the managed
PostgreSQL layer during backend refresh, with explicit constraints for
`AiiDA/Postgres` ownership boundaries.

## Decision Summary

Recommendation: **Adopt Neon with staged rollout guardrails**.

- Adopt now for non-execution-critical durable stores (projection-side metadata,
  adapter-side operational metadata).
- Adopt for AiiDA metadata only after staging validation gates pass.
- Do not block backend refresh on immediate full Neon migration.

This keeps momentum while preventing hidden durability regressions in
execution-critical paths.

## Connection and Migration Constraints

1. Endpoint selection
   - Use Neon direct connection endpoints for schema migration operations and
     long-running maintenance.
   - Use pooled endpoints for short-lived API request traffic where transaction
     pooling is acceptable.
2. Session semantics
   - Avoid relying on session-level connection state in pooled mode.
   - Run migration and admin commands through direct endpoints only.
3. Autosuspend and cold-start handling
   - Account for compute resume latency on idle databases.
   - Keep execution-critical environments on non-aggressive suspend settings
     during active windows.
4. Connection budgeting
   - Set application-level pool caps per service (FastAPI adapter, workers,
     background jobs) to stay below plan limits with headroom.
   - Treat Neon plan connection/storage limits as hard SLO boundaries.

## AiiDA/Postgres Compatibility Assumptions

1. Backend model
   - AiiDA production usage expects PostgreSQL-backed storage (`core.psql_dos`)
     and a message broker.
   - SQLite defaults are for local experimentation, not production durability.
2. Storage split
   - AiiDA durability is not database-only: metadata in PostgreSQL plus
     repository/object-store contents must stay consistent.
3. Migration discipline
   - Keep AiiDA schema migrations serialized and rehearsed against a Neon branch
     clone before applying to higher environments.
4. Operational ownership
   - Continue ADR-0001 ownership split:
     - execution internals -> AiiDA/Postgres
     - product projections -> Convex
   - Neon is a PostgreSQL substrate choice, not an ownership change.

## Operational Guardrails

1. Backup and restore
   - Use Neon restore/branch recovery for database rollback.
   - Pair every database backup policy with repository/object-store backup policy
     so AiiDA metadata and file artifacts can be recovered together.
   - Run periodic restore drills and record RTO/RPO outcomes.
2. Environment strategy
   - One long-lived branch per environment (`dev`, `staging`, `prod`).
   - Short-lived feature branches allowed only for migration rehearsal and
     destructive test scenarios.
3. Capacity and cost controls
   - Monitor Neon compute usage, storage growth, and active connections.
   - Define alert thresholds before saturation (for example 70/85/95% bands).
4. Release safety
   - For each migration PR: include forward migration, rollback path,
     and data consistency checks.
   - No dual-write fallback as steady state; use short-lived operational rollback
     only.

## Adopt/Defer Recommendation

Decision: **Adopt (staged)**.

- Adopt in current cycle for staging and non-critical production data paths.
- Gate execution-critical AiiDA production cutover on the following criteria:
  1. connection budget test passes under expected daemon/API concurrency
  2. migration rehearsal succeeds on a production-like Neon branch clone
  3. combined DB + repository restore drill meets target RTO/RPO
  4. one full release cycle without Neon-specific availability incidents in
     staging

If any gate fails, defer only the AiiDA production cutover while continuing
backend refresh on other tracks.

## Source Notes (validated 2026-02-15)

- Neon docs:
  - <https://neon.com/docs/connect/connection-pooling>
  - <https://neon.com/docs/introduction/branching>
  - <https://neon.com/docs/introduction/plans>
  - <https://neon.com/docs/guides/backup-restore>
- AiiDA docs:
  - <https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_quick.html>
  - <https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/installation.html>
