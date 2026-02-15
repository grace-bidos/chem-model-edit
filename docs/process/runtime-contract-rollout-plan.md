# Runtime Contract Rollout Plan

This document turns ADR-0004 and runtime contracts into merge-ready parallel delivery lanes.

## Lane split

- Lane A (`Show`, merge step 4): Convex projection enrichment and monotonicity enforcement
- Lane B (`Show`, merge step 3): FastAPI command/event/projection contract enforcement
- Lane C (`Show`, merge step 2): Management-node integration and event emission compatibility
- Lane D (`Ship`, merge step 5): Redis/RQ retirement cleanup and guardrail tests

## Merge order

1. Contract-first docs PR (this phase)
2. Lane C: management-node integration and event emission compatibility
3. Lane B: FastAPI contract enforcement (`SubmitJobCommand` and `ExecutionEvent` validation)
4. Lane A: Convex projection monotonicity enforcement
5. Lane D: Redis/RQ retirement path removals and guardrail tests

## Required checks per lane

- Contract docs referenced in PR description
- SoR mapping table for changed fields
- Idempotency test coverage for changed boundary
- Tenant boundary tests for changed boundary

## Exit criteria for contract rollout milestone

- `ADR-0004` accepted and linked in implementation PRs
- `docs/contracts/*` referenced by all runtime lane PRs
- No open PR with ambiguous ownership across Convex, FastAPI, and management node
