# Contract: Slurm Real-Adapter Cutover Split

## Purpose

Define the contract for cutting over from the current internal Slurm
stub boundary (`slurm-adapter-boundary/v1`) to a real Slurm adapter
implementation without changing Product-facing runtime contracts.

This contract applies ADR-0003 single-write cutover policy to the Slurm adapter
lane: no long-lived dual-write semantics, rollback is operational and
short-lived, and ownership boundaries remain explicit.

## Scope

- In scope: internal adapter boundary behavior in FastAPI runtime path.
- In scope: cutover controls, preconditions, and acceptance gates.
- Out of scope: Slurm infrastructure setup or provisioning execution.
- Out of scope: changes to `SubmitJobCommand`, `ExecutionEvent`, or
  `ProjectionUpdate` payload schemas.

## Producer and consumer

- Producer: `FastAPI` orchestration adapter (`services/zpe/slurm_policy.py` path)
- Consumer: submission runtime path in `FastAPI` (`services/zpe/backends.py`)
- Product/Execution contract context:
  - Product SoR: `Convex`
  - Execution SoR: `management node stack (AiiDA/Postgres + Slurm)`

## Contract baseline (current)

Current internal boundary emits `SlurmAdapterBoundaryStub` with:

- `adapter`: `"passthrough"` or `"stub-policy"`
- `contract_version`: `slurm-adapter-boundary/v1`
- `requested_queue`, `resolved_queue`, `used_fallback`
- `mapping` (nullable) with:
  - `queue`, `partition`, `account`, `qos`, `max_walltime_minutes`

## Interface delta contract for real-adapter cutover

Required deltas for implementation lanes:

- Keep `contract_version` as `slurm-adapter-boundary/v1` during cutover.
- Add adapter mode `"real-policy"` behind the same boundary shape.
- Preserve existing queue resolution fields and `mapping` keys exactly.
- Keep `SubmitJobCommand`, `ExecutionEvent`, and `ProjectionUpdate` schemas
  unchanged.
- Allow additive internal metadata only if it is optional and does not alter
  current retry/idempotency behavior.

Non-allowed deltas in this cutover:

- Changing existing field names/types in boundary response.
- Coupling cutover to external API/OpenAPI schema changes.
- Introducing long-lived parallel ownership between stub and real adapters.

## Idempotency key

No idempotency key shape changes are allowed in this cutover. Existing
submission idempotency contract remains authoritative.

## Retry policy

Retry behavior must remain equivalent across `stub-policy` and `real-policy`.
Any new adapter-local retry loop must be bounded and must not change existing
runtime retry semantics.

## Error semantics

- Invalid or missing adapter policy/config remains a configuration error.
- Queue denial remains an explicit policy denial outcome.
- Real-adapter backend unavailability must produce operational failure signals
  without silently mutating queue resolution semantics.

## Infrastructure preconditions

Before enabling `"real-policy"` in any runtime environment:

- Management node has reachable Slurm control plane for runtime submission path.
- Adapter runtime credentials/config are present and readable.
- Policy source is valid and compatible with queue mapping requirements.
- Health/readiness checks prove adapter dependency availability.
- Runbook for rollback switch is published and linked in PR description.

## Rollback switches (operational, short-lived)

Primary rollback controls:

- `ZPE_SLURM_ADAPTER=stub-policy` (fallback to policy-backed stub path)
- `ZPE_SLURM_ADAPTER=passthrough` (emergency-safe direct queue passthrough)

Rollback rules:

- Rollback flags are emergency-stop controls with explicit owner and expiry.
- Rollback does not change SoR boundaries (`Convex` Product SoR,
  management-node Execution SoR).
- After incident stabilization, rollback settings must be removed or reset in
  a tracked follow-up task.

## Tenant and security boundary

- Adapter mode changes must not alter tenant scoping rules.
- No cross-tenant queue resolution state may be introduced.
- Credentials used for real-adapter access must be environment-scoped and
  non-user-authoritative.

## Observability fields

At minimum, cutover and rollback telemetry should include:

- `trace_id`
- `tenant_id`
- `workspace_id` (when available)
- `job_id` (when available)
- `slurm_adapter`
- `slurm_contract_version`
- `requested_queue_name`
- `resolved_queue_name`
- `used_fallback`

## SoR mapping

- Product-visible state owner remains `Convex` (Product SoR).
- Execution internals owner remains management node stack
  (`AiiDA/Postgres + Slurm`) (Execution SoR).
- FastAPI remains orchestration adapter and does not become a durable SoR.

## Acceptance criteria

- AC-001: Real adapter can be enabled via config switch without changing any
  public runtime contract schema.
- AC-002: Existing stub-path tests remain green for `passthrough` and
  `stub-policy`.
- AC-003: New tests verify `real-policy` keeps boundary field semantics
  (`requested_queue`, `resolved_queue`, fallback semantics, mapping fields).
- AC-004: Contract docs and implementation stay in lockstep in same PR layer.
- AC-005: Rollback procedure is tested as a config-only operation.
- AC-006: No PR in this cutover lane introduces Redis/RQ business ownership.
- AC-007: Tenant/security behavior remains unchanged across adapter modes.

## Follow-up implementation split proposal (child-task checklist)

- [ ] Implement real Slurm adapter mode behind existing boundary (`Show`, size M)
      Title: `Implement real-policy Slurm adapter mode behind v1 boundary`
- [ ] Add startup/readiness precondition enforcement for real adapter
      (`Show`, size S)
      Title: `Enforce real-adapter preconditions in runtime readiness`
- [ ] Add cutover and rollback observability fields/logging (`Ship`, size S)
      Title: `Add real-adapter cutover and rollback observability contract checks`
- [ ] Expand contract/integration tests for stub + real parity (`Show`, size M)
      Title: `Validate stub-real adapter parity and rollback contract`
- [ ] Remove temporary rollback exception metadata after stabilization
      (`Ship`, size XS)
      Title: `Retire temporary rollback switches for Slurm adapter cutover`

## Verification checklist for PR review

- Does the PR explicitly state adapter mode transition (`stub-policy` ->
  `real-policy`)?
- Are boundary field semantics unchanged and documented?
- Are rollback controls and expiry/removal ownership documented?
- Are tenant/security and idempotency/retry semantics unchanged?
- Is this contract doc updated in the same PR as behavior changes?
