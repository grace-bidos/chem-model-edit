# ADR-0002: Canonical JobState Contract and Transition Rules

- Status: Accepted
- Date: 2026-02-15
- Linear: `GRA-13`
- Owners: Backend platform

## Context

Job state values were historically stringly-typed and mutable from several code paths. Backend refresh needs a canonical, shared, and enforceable state machine to avoid invalid transitions during migration.

## Decision

Canonical state set:

- `queued`
- `started`
- `finished` (terminal)
- `failed` (terminal, with explicit retry path to `queued`)

Allowed transitions:

- `<none> -> queued|started|finished|failed`
- `queued -> started|failed|finished|queued`
- `started -> queued|failed|finished|started`
- `failed -> queued|failed`
- `finished -> finished`

Disallowed transitions:

- Any transition from `finished` to non-`finished`.
- `failed -> started|finished` without first re-queueing.

## Consequences

- State values are validated centrally before persistence.
- API schema for job status now exposes the canonical literal union, not free-form strings.
- Invalid transitions fail fast with explicit errors.

## Implementation

Implemented in this change:

- `apps/api/services/zpe/job_state.py`: canonical `JobState` type + transition validation.
- `apps/api/services/zpe/result_store.py`: transition guard in `set_status`.
- `apps/api/app/schemas/zpe.py`: `ZPEJobStatus.status` uses canonical `JobState`.
- `packages/shared/src/job-state.ts`: shared TS contract for follow-up web/BFF integration.

## Verification

- Unit tests ensure:
  - invalid transition rejection from `finished` to `started`
  - valid retry transition from `failed` to `queued`
- Existing compute result paths continue to pass and respect terminal-state protection.
