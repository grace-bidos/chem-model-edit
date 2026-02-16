# AiiDA Promotion Gate: sqlite Smoke -> PostgreSQL + RabbitMQ VM Runtime (GRA-125)

This Ask deliverable defines a promotion gate to move from local sqlite smoke confidence to VM-backed `core.psql_dos` runtime confidence before enabling PostgreSQL + RabbitMQ AiiDA runtime as the active execution path.

## Scope and decision

In scope:

- Gate the runtime promotion decision only.
- Reuse existing bootstrap/smoke scripts.
- Require evidence that is reproducible from commands and logs.

Out of scope:

- Runtime behavior changes.
- Contract/schema changes.

Promotion decision outputs:

- `PROMOTE`: VM runtime is eligible for the next Show lane.
- `HOLD`: stay on sqlite smoke baseline and execute fallback actions.

## Gate inputs

Required prerequisites:

- VM bootstrap baseline completed via `docs/process/aiida-vm-bootstrap.md`.
- Slurm connectivity smoke path available via `docs/process/aiida-slurm-smoke-vm.md`.
- Runtime contract boundary remains aligned with `docs/adr/ADR-0001-system-of-record-boundaries.md` and `docs/adr/ADR-0004-runtime-contract-boundary-byo-aiida.md`.

## Executable checklist

Run in a clean worktree on the target branch.

### 1) Environment and tool preflight

```bash
for tool in node corepack pnpm uv; do command -v "$tool" >/dev/null || { echo "missing: $tool"; exit 1; }; done
node -v
pnpm -v
uv --version
```

Pass threshold:

- all commands return exit `0`

Fallback criteria:

- any missing tool or non-zero exit -> `HOLD`

Required evidence:

- command transcript captured in `investigations/artifacts/gra-125/<timestamp>/01-preflight.log`

### 2) sqlite smoke baseline (control)

Run current sqlite smoke path exactly as-is for comparison.

```bash
# Replace with the project-standard sqlite smoke command used by lane owner.
# Keep stdout/stderr.
```

Pass threshold:

- 3 consecutive runs succeed (`exit 0`) with no contract/test failures
- p95 wall time does not regress by more than 25% versus previous known-good sqlite smoke median

Fallback criteria:

- fewer than 3/3 successes, or p95 regression > 25% -> `HOLD`

Required evidence:

- `02-sqlite-smoke-run1.log`
- `02-sqlite-smoke-run2.log`
- `02-sqlite-smoke-run3.log`
- `02-sqlite-smoke-summary.md` (success count + timing table)

### 3) VM bootstrap verification (PostgreSQL + RabbitMQ + AiiDA profile)

```bash
scripts/aiida-vm-bootstrap.sh \
  --env-file ops/aiida-vm/aiida-vm.env \
  --sanity-check
```

If first run on the VM for this profile:

```bash
scripts/aiida-vm-bootstrap.sh \
  --apply \
  --env-file ops/aiida-vm/aiida-vm.env \
  --init-profile \
  --sanity-check
```

Pass threshold:

- sanity-check exits `0`
- `verdi profile list` contains target profile
- `verdi status` shows daemon and broker reachable

Fallback criteria:

- missing profile, broker unreachable, or sanity-check non-zero -> `HOLD`

Required evidence:

- `03-vm-bootstrap.log`
- `03-verdi-profile-list.log`
- `03-verdi-status.log`

### 4) AiiDA -> Slurm VM smoke

```bash
scripts/aiida-slurm-smoke-vm.sh --artifact-dir investigations/artifacts/gra-125/<timestamp>/slurm-smoke
```

Pass threshold:

- script exits `0`
- `scontrol ping` reports controller `UP`
- `sinfo` shows at least one schedulable node
- `verdi computer test` passes for target `core.slurm` computer

Fallback criteria:

- script exits `2` (blocked) or any required check missing -> `HOLD`

Required evidence:

- full artifact directory from the script run
- `04-slurm-smoke-summary.md` with pass/fail mapping

### 5) FastAPI/Convex validation-flow compatibility gate

This promotion gate must not blur ownership boundaries.

Run boundary/contract checks:

```bash
PYTHONPATH=apps/api uv run --project apps/api python apps/api/scripts/export_openapi.py --check
pnpm -C packages/api-client run generate
git diff --exit-code -- packages/api-client/src/generated/schema.ts
uv run pytest apps/api/tests/test_convex_event_relay.py
```

Pass threshold:

- no OpenAPI/client artifact drift
- Convex projection relay test suite passes
- no new field ownership ambiguity against ADR-0001/ADR-0004 principles

Fallback criteria:

- contract drift or failing Convex relay tests -> `HOLD`

Required evidence:

- `05-openapi-check.log`
- `05-api-client-generate.log`
- `05-convex-relay-tests.log`
- `05-boundary-notes.md` documenting ownership check outcome

### 6) Promotion decision record

Create final decision artifact:

- `investigations/artifacts/gra-125/<timestamp>/promotion-decision.md`

Include:

- decision: `PROMOTE` or `HOLD`
- gate checklist result table (1-5)
- known risks and mitigations
- operator + reviewer
- date/time (UTC)

Pass threshold:

- all gates 1-5 pass and evidence bundle is complete

Fallback criteria:

- any gate failure or missing evidence file -> `HOLD`

## Minimum acceptance thresholds (summary)

- Stability: sqlite smoke 3/3 success before promotion call.
- Runtime readiness: VM bootstrap sanity-check and AiiDA/Slurm smoke both pass.
- Contract safety: OpenAPI/client drift checks pass and Convex relay boundary tests pass.
- Evidence quality: every gate has machine-generated logs plus short operator summary.

## Fallback and rollback policy

When gate result is `HOLD`:

- Keep sqlite smoke as release-blocking confidence source.
- Do not enable PostgreSQL+RabbitMQ runtime as default execution path.
- Open a follow-up child issue with:
  - failing gate step
  - copied evidence paths
  - explicit unblock condition
- Re-run only the failed gates after remediation; retain previous passing evidence.

## Implications for FastAPI/Convex validation flow

- FastAPI remains orchestration/validation adapter, not durable owner.
- Convex remains product-facing projection SoR; promotion does not permit broad AiiDA-internal mirroring.
- AiiDA/PostgreSQL remains execution-internal SoR for runtime/provenance details.
- Promotion is blocked if checks imply contract drift or ownership drift across `SubmitJobCommand`, `ExecutionEvent`, or `ProjectionUpdate` boundaries.

Reference contracts:

- `docs/contracts/command-submit-job.md`
- `docs/contracts/event-execution-lifecycle.md`
- `docs/contracts/projection-update.md`
