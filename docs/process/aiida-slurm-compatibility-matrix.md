# AiiDA/Slurm Compatibility Matrix and Upgrade Canary Policy (GRA-132)

This document defines the executable compatibility policy for the VM runtime path:

- Python runtime in `apps/api`
- `aiida-core`
- Slurm control plane
- PostgreSQL
- RabbitMQ

Use this as the source of truth for promotion-gate canary decisions before runtime cutover.

## Scope and authority

- Scope: VM runtime validation and promotion gating (`docs/process/aiida-promotion-gate-vm.md`).
- Gate mapping: Tier-2 and Tier-3 from `docs/process/slurm-verification-boundary.md`.
- Applies to both local VM runs and remote VM reruns (`scripts/aiida-promotion-gate-vm-remote.sh`).

## Supported ranges (project policy)

This table defines the currently supported window for gate decisions.

| Component | Supported range | Source/basis | Gate level required |
| --- | --- | --- | --- |
| Python | `3.13.x` | `apps/api/pyproject.toml` requires `>=3.13`; runtime images are `python:3.13-slim` | Tier-2 required |
| AiiDA (`aiida-core`) | `2.7.x` | dependency group in `apps/api/pyproject.toml` (`>=2.7.3`) plus current VM scripts | Tier-2 required |
| Slurm | `23.11.x` to `24.05.x` | current VM smoke/tooling assumptions (`scontrol`, `sinfo`, `core.slurm`) | Tier-2 required |
| PostgreSQL | `14.x` to `16.x` | bootstrap flow installs distro package and uses `core.psql_dos` profile | Tier-2 required |
| RabbitMQ | `3.12.x` to `3.13.x` | bootstrap flow uses `rabbitmqctl` user/vhost/permission contract | Tier-2 required |

Policy outside supported range:

- Same-major patch/minor bump outside table but without known breakage: `Tier-3 canary required` before promotion.
- Major upgrade (`x -> x+1`): `Tier-2 + Tier-3 canary both required`; default decision is `HOLD` until both pass.

## Known unsafe combinations (block promotion)

| Combination | Unsafe reason | Typical signal | Required action |
| --- | --- | --- | --- |
| Python `3.13.x` + AiiDA `<2.7` | outside tested project contract | `uv sync`/`verdi` setup failures | pin AiiDA back to `2.7.x` or complete full validation lane |
| RabbitMQ `4.x` + current AiiDA VM path | unvalidated major bump for broker behavior in this project | `verdi status` broker/daemon errors or bootstrap permission drift | rollback RabbitMQ to `3.12/3.13` |
| PostgreSQL `17.x` + current VM profile path | unvalidated major bump for project runtime gate | profile setup/migration anomalies during gate 3 | rollback PostgreSQL to `14-16` |
| Slurm `25.x` + current smoke scripts | unvalidated scheduler major bump | `verdi computer test` scheduler plugin failures | rollback Slurm to `23.11/24.05` or run dedicated validation lane |
| Mixed unsupported majors (2 or more components) | root-cause isolation impossible in canary | multiple failing gate artifacts in one run | rollback to previous known-good bundle first, then re-test one component at a time |

## Version capture (mandatory before any upgrade canary)

Run on target VM and store under canary artifacts:

```bash
mkdir -p investigations/artifacts/gra-132/<timestamp>/versions
python3 --version | tee investigations/artifacts/gra-132/<timestamp>/versions/python.log
uv run --project apps/api --group aiida verdi --version | tee investigations/artifacts/gra-132/<timestamp>/versions/verdi.log
scontrol -V | tee investigations/artifacts/gra-132/<timestamp>/versions/slurm.log
psql --version | tee investigations/artifacts/gra-132/<timestamp>/versions/postgres.log
rabbitmqctl version | tee investigations/artifacts/gra-132/<timestamp>/versions/rabbitmq.log
```

## Upgrade canary procedure (executable)

### Phase A: baseline control (current known-good bundle)

1. Run gate 3 and gate 4 using existing runbook:
   - `scripts/aiida-vm-bootstrap.sh --env-file ops/aiida-vm/aiida-vm.env --sanity-check`
   - `scripts/aiida-slurm-smoke-vm.sh --artifact-dir investigations/artifacts/gra-132/<timestamp>/baseline`
2. Result must be `exit 0` for both; otherwise upgrade canary is not allowed to start.

### Phase B: single-component canary

1. Change exactly one component version (for example PostgreSQL `15 -> 16`).
2. Re-run:
   - bootstrap sanity (`gate 3`)
   - Slurm smoke (`gate 4`)
3. Keep all artifacts under:
   - `investigations/artifacts/gra-132/<timestamp>/<component>-canary/`

### Phase C: decision

- `PROMOTE_CANDIDATE` when all are true:
  - gate 3 `exit 0`
  - gate 4 `exit 0`
  - no blocker diagnostics generated
  - no contract drift signals from promotion gate step 5
- otherwise `ROLLBACK_REQUIRED`.

## Rollback criteria (hard stop)

Rollback immediately to previous known-good versions when any condition is met:

- `scripts/aiida-vm-bootstrap.sh --sanity-check` exits non-zero.
- `scripts/aiida-slurm-smoke-vm.sh` exits `2` (blocked) after upgrade.
- `verdi status` shows broker/daemon not reachable for upgraded stack.
- `verdi -p <profile> computer test <label>` fails after previously passing baseline.
- Promotion-gate step 5 reports OpenAPI/client drift or Convex relay compatibility failure caused by runtime dependency changes.

Rollback completion checklist:

1. Restore previous package/runtime versions.
2. Re-run baseline gate 3 and gate 4.
3. Confirm both return `exit 0`.
4. Attach rollback evidence and root-cause note to the issue.

## Tier-2 / Tier-3 gate mapping

| Change type | Tier-2 requirement | Tier-3 requirement | Promotion rule |
| --- | --- | --- | --- |
| Patch/minor within supported range | required | optional | Tier-2 pass -> may promote |
| Minor outside supported range but same major | required | required | promote only after both pass |
| Any major upgrade | required | required | default `HOLD` until both pass and rollback plan rehearsed |
| Multiple components changed together | required (after each isolated canary) | required | only promote after isolated canaries + combined run pass |

## Operator notes for reproducibility

- Change one component per canary run. Do not batch multiple major upgrades.
- Reuse existing artifact naming conventions (`gra-125`, `gra-126`, `gra-132`) and keep UTC timestamps.
- For remote execution, use `scripts/aiida-promotion-gate-vm-remote.sh` to collect deterministic evidence for gates 3/4.
