# Contract: ExecutionEvent

## Purpose

Defines lifecycle events sent from management node (`AiiDA + Slurm`) to FastAPI for state progression.

## Producer and consumer

- Producer: management node event emitter
- Consumer: FastAPI orchestration adapter

## Required fields

- `event_id` (string, globally unique)
- `tenant_id` (string)
- `workspace_id` (string)
- `job_id` (string)
- `submission_id` (string)
- `execution_id` (string, execution SoR id; stable id for the AiiDA/Slurm execution instance)
- `state` (enum)
  - `accepted` (execution request acknowledged by execution owner)
  - `running` (runtime started and progressing)
  - `completed` (successful terminal state)
  - `failed` (failed terminal state)
- `occurred_at` (RFC3339 timestamp)
- `trace_id` (string)

## Optional fields

- `status_detail` (string)
- `scheduler_ref` (object)
  - `slurm_job_id` (string)
  - `partition` (string)
  - `qos` (string)
- `result_ref` (object, required when `state=completed`)
  - `output_uri` (string)
  - `metadata_uri` (string, optional)
- `error` (object, required when `state=failed`)
  - `code` (string)
  - `message` (string)
  - `retryable` (boolean)

## Idempotency key

- `event_id` is the event idempotency key.
- Duplicate `event_id` must be treated as idempotent replay and return success with no additional side effects.

## Retry policy

- Producer retries delivery on timeout or non-2xx.
- Consumer must be at-least-once safe using `event_id` dedup.

## State monotonicity

- Allowed transitions:
  - `accepted -> running -> completed`
  - `accepted -> running -> failed`
  - `accepted -> failed`
- Non-monotonic regressions are rejected as `409 invalid_transition`.
- State-specific required fields:
  - `state=completed`: `result_ref.output_uri` is required.
  - `state=failed`: `error.code`, `error.message`, and `error.retryable` are required.

## Error semantics

- `400 invalid_event`: schema/required-field failure
- `403 forbidden`: tenant boundary violation
- `404 submission_not_found`: unknown submission reference
- `409 invalid_transition`: non-monotonic or conflicting state

## Tenant and security boundary

- `tenant_id` must match submission owner before applying event.
- Events are scoped to one tenant; cross-tenant event references are rejected.

## Observability fields

- Required in logs: `trace_id`, `tenant_id`, `job_id`, `submission_id`, `execution_id`, `event_id`, `state`.

## SoR mapping

- Execution lifecycle facts (`execution_id`, `state`, `scheduler_ref`, `result_ref`) are Execution SoR (`management node stack`, `AiiDA/Postgres + Slurm`).
- Product-facing references (`tenant_id`, `workspace_id`, `job_id`, `submission_id`) remain Product SoR joins.
- FastAPI consumes and validates events but does not own execution durability.

## Compatibility notes (AiiDA/Slurm identifiers)

- `execution_id` should be a stable execution-owner identifier when available from AiiDA (for example, process/workchain id).
- During migration, when a dedicated AiiDA execution id is not yet emitted, producers may use a deterministic fallback scoped to the management node and lease (for example `management_node_id:job_id:lease_id`) and keep it stable for replay of the same lifecycle event.
- `scheduler_ref.slurm_job_id` is optional because Slurm allocation can occur after earlier lifecycle states; when present it should remain immutable for the same execution.
- `scheduler_ref.partition` and `scheduler_ref.qos` map to Slurm scheduling metadata and are optional enrichment fields.
