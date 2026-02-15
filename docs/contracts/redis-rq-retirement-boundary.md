# Contract: Redis/RQ Retirement Boundary

## Purpose

Defines what Redis/RQ may still do during transition and what is forbidden, so retirement work can proceed without semantic drift.

## Producer and consumer

- Transitional producers: FastAPI adapter and remaining legacy worker paths during migration only
- Transitional consumers: migration-safe runtime adapters that use bounded TTL cache/lease keys
- Final target consumers/producers: runtime contract surfaces in `docs/contracts/` without RQ business dependency

## Required fields

This retirement boundary is policy-oriented and does not define a runtime payload schema. Required governance fields for exceptions are:

- `owner` (team/person)
- `exception_reason`
- `expiry_date`
- `removal_issue`

## Idempotency key

Not a payload contract. Idempotency for runtime behavior is defined in:

- `command-submit-job.md`
- `event-execution-lifecycle.md`
- `projection-update.md`

## Retry policy

- Runtime retries must not depend on Redis/RQ owning business truth.
- Retries remain valid only when durable owners (`Convex` or `management node stack`) can safely replay.

## Allowed transitional usage

- Redis as short-lived cache with explicit TTL
- Redis as transient lease/lock transport
- Redis as ephemeral dedup buffer for in-flight processing only

## Forbidden usage

- Redis as durable owner for business fields (`job status`, `result metadata`, `tenant ownership`, `queue policy`)
- New runtime features built on RQ execution pipelines
- Dual-write designs that keep Redis and durable owners in long-lived semantic parity

## Migration cut lines

- Product-visible state owner: Convex only
- Execution internals owner: management node (`AiiDA/Postgres + Slurm`) only
- FastAPI role: routing/validation/orchestration only

## Error semantics

- PR-level policy violation should block merge if Redis/RQ business ownership is introduced or expanded without approved temporary exception metadata.

## Tenant and security boundary

- Transitional Redis usage must remain tenant-scoped and non-authoritative.
- No cross-tenant business authority may be derived from Redis keys.

## Observability fields

For remaining transitional Redis keys/logs, include at minimum:

- `trace_id`
- `tenant_id`
- `workspace_id` (when applicable)
- `job_id` (when applicable)
- `key_ttl_seconds`

## SoR mapping

- Redis/RQ are non-SoR transitional infrastructure only.
- Product SoR: Convex.
- Execution SoR: management node stack (`AiiDA/Postgres + Slurm`).

## Acceptance criteria for retirement completion

- No write path where RQ workers determine final business state
- No user-visible endpoint reading authoritative job state from Redis
- All runtime paths covered by contract surfaces in `docs/contracts/`
- Redis keys used only for bounded-TTL operational concerns

## Enforcement checklist for PR review

- Did this PR add or expand Redis business semantics?
- Did this PR introduce new RQ runtime dependency?
- Are ownership and contract docs updated for changed runtime behavior?

If any answer is yes to the first two, PR must be rejected unless this document is explicitly updated with a temporary expiry-bound exception.
