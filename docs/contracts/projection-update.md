# Contract: ProjectionUpdate

## Purpose

Defines the update contract emitted by FastAPI to Convex for product-facing projection state.

## Producer and consumer

- Producer: FastAPI orchestration adapter
- Consumer: Convex projection mutation endpoint

## Required fields

- `projection_event_id` (string)
- `tenant_id` (string)
- `workspace_id` (string)
- `job_id` (string)
- `submission_id` (string)
- `projection_state` (enum)
  - `queued`
  - `running`
  - `succeeded`
  - `failed`
- `source_execution_state` (enum)
  - `accepted`
  - `running`
  - `completed`
  - `failed`
- `occurred_at` (RFC3339 timestamp)
- `trace_id` (string)

## Optional fields

- `display_message` (string)
- `result_ref` (object)
  - `output_uri` (string)
  - `summary` (string, optional)
- `failure` (object)
  - `code` (string)
  - `message` (string)

## Idempotency key

- `projection_event_id` is dedup key at Convex boundary.
- Duplicate updates with same key must be treated as no-op.

## Retry policy

- FastAPI retries projection writes on transient failures.
- Retries must preserve `projection_event_id` to avoid duplicate user-visible updates.

## Projection monotonicity

- Monotonic order for user-visible state:
  - `queued -> running -> succeeded|failed`
- Regressive updates with a new `projection_event_id` must be rejected with `409 invalid_projection_transition`.
- Exact duplicate `projection_event_id` replays must be treated as idempotent no-op success.

## Error semantics

- `400 invalid_projection_update`: schema violation
- `403 forbidden`: tenant mismatch
- `404 projection_target_not_found`: job projection missing
- `409 invalid_projection_transition`: non-monotonic update

## Tenant and security boundary

- Projection writes are always tenant-scoped.
- Consumer must validate tenant ownership before mutation.

## Observability fields

- Required: `trace_id`, `tenant_id`, `workspace_id`, `job_id`, `submission_id`, `projection_event_id`, `projection_state`.

## SoR mapping

- User-visible projection state (`projection_state`, projection payload) is Product SoR (`Convex`).
- Execution-origin references (`source_execution_state`, `result_ref`, `failure`) are derived from Execution SoR events.
- FastAPI is the adapter that translates execution events into projection updates.
