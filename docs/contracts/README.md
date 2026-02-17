# Runtime Contracts

This directory defines the fixed runtime contracts for BYO AiiDA/Slurm execution.

## Contract set

- `command-submit-job.md`: command contract for job submission to FastAPI.
- `event-execution-lifecycle.md`: event contract from management node to FastAPI.
- `projection-update.md`: projection update contract from FastAPI to Convex.
- `redis-rq-retirement-boundary.md`: migration boundary and retirement gate for Redis/RQ.
- `slurm-real-adapter-cutover.md`: stub-to-real Slurm adapter cutover contract and acceptance gates.

## Contract policy

Every runtime contract document must include:

- Purpose
- Producer and consumer
- Required fields
- Idempotency key
- Retry policy
- Error semantics
- Tenant/security boundary
- Observability fields
- SoR mapping

Breaking changes are allowed in current pre-production phase, but contract docs must be updated in the same PR as code changes.
