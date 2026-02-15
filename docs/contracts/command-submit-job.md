# Contract: SubmitJobCommand

## Purpose

Defines the command accepted by FastAPI to submit an execution request for a registered management node.

## Producer and consumer

- Producer: BFF/Convex-facing backend caller (authenticated tenant context)
- Consumer: FastAPI orchestration adapter

## Required fields

- `tenant_id` (string)
- `workspace_id` (string)
- `job_id` (string, product-facing stable id)
- `idempotency_key` (string)
- `management_node_id` (string)
- `execution_profile` (object)
  - `queue_name` (string)
- `resource_shape` (object)
  - `cpu` (integer, vCPU count)
  - `memory_mib` (integer, MiB)
  - `walltime_seconds` (integer, seconds)
- `payload_ref` (object)
  - `input_uri` (string, URI such as `s3://...`, `gs://...`, or signed `https://...`)
- `requested_by` (object)
  - `user_id` (string)

## Optional fields

- `execution_profile.qos` (string)
- `execution_profile.account` (string)
- `payload_ref.artifact_bucket` (string, storage location name or URI)
- `requested_by.session_id` (string)

## Idempotency key

- `idempotency_key` must uniquely identify a submission intent within (`tenant_id`, `workspace_id`).
- Replays with same key return previous accepted response only when all required fields match exactly:
  - `tenant_id`, `workspace_id`, `job_id`, `management_node_id`, `execution_profile`, `payload_ref.input_uri`, `requested_by.user_id`.
- Optional fields may differ only for client hints and must not change execution intent:
  - `requested_by.session_id`.
- Any differences in other optional fields are treated as intent changes and rejected with conflict.
- Replays with same key and different payload return `409 idempotency_conflict`.

## Retry policy

- Client may retry on transport errors and `5xx`.
- FastAPI must treat retries as idempotent using `idempotency_key`.

## Success response

- `202 Accepted`
- Body:
  - `job_id`
  - `submission_id` (server-generated unique opaque identifier, immutable)
  - `execution_owner` (`management_node_id`)
  - `accepted_at`
  - `trace_id`
- `submission_id` constraints:
  - Treat as opaque string.
  - Recommended format: UUIDv4-compatible identifier.
  - Used for status query correlation and audit linkage.

## Error semantics

- `400 invalid_request`: schema/validation failure
- `401 unauthorized`: auth/session invalid
- `403 forbidden`: tenant/workspace access violation
- `404 management_node_not_found`: unknown or inaccessible node
- `409 idempotency_conflict`: reused key with different intent
- `422 policy_rejected`: queue/account/qos/resource policy mismatch
- `503 execution_unavailable`: management node temporarily unavailable

## Tenant and security boundary

- `tenant_id` is mandatory and validated before any downstream call.
- `management_node_id` must be tenant-owned or explicitly shared to tenant.
- No cross-tenant read or write side effects are allowed.

## Observability fields

- Required in logs/events: `trace_id`, `tenant_id`, `workspace_id`, `job_id`, `submission_id`, `management_node_id`, `idempotency_key`.

## SoR mapping

- Product-facing identifiers (`tenant_id`, `workspace_id`, `job_id`) are Product SoR (`Convex`).
- Execution target and runtime details are Execution SoR (`management node stack`, `AiiDA/Postgres + Slurm`).
- FastAPI does not become durable owner; it validates and routes command processing.
