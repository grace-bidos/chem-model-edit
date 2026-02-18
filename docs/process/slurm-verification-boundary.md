# Slurm Verification Boundary (GRA-66)

This note defines what we should verify without a real Slurm cluster vs what must run against real Slurm.

## Tier Definitions

- Tier-1 (contract/mocks): deterministic checks runnable on developer machines/CI without `slurmctld`.
- Tier-2 (local VM real Slurm smoke): minimal real-runtime checks on a VM where Slurm + AiiDA are installed.
- Tier-3 (optional shared cluster): non-blocking confidence checks on a shared/remote cluster environment.

## Boundary Map

| Check area | Tier | Why this tier | Concrete repo command examples |
| --- | --- | --- | --- |
| Slurm queue policy contract and fallback semantics (`route-default` / `deny`) | Tier-1 | Pure Python/JSON validation; no scheduler RPC | `uv run pytest apps/api/tests/test_zpe_slurm_policy.py` |
| Real-adapter precondition readiness contract (`real-policy`, rollback guard, policy-file validation) | Tier-1 | FastAPI readiness semantics and adapter-mode gating can be validated without running full scheduler stack | `uv run pytest apps/api/tests/test_health_api.py -k real_adapter` |
| BYO onboarding schema and queue-resolution CLI contract | Tier-1 | Script + fixtures only; no real `sinfo`/`sbatch` invocation | `uv run pytest apps/api/tests/test_validate_slurm_onboarding.py`<br>`python3 scripts/validate-slurm-onboarding.py --policy apps/api/config/slurm/onboarding.policy.example.json --registration apps/api/config/slurm/node-registration.example.json --requested-queue standard` |
| Runtime event payload compatibility carrying Slurm metadata fields (`slurm_job_id`, partition, qos) | Tier-1 | Runtime contract validation without scheduler behavior dependency | `uv run pytest apps/api/tests/test_runtime_api.py` |
| Queue target ownership and selection behavior (`fakeredis`) | Tier-1 | Redis mocked; no Slurm runtime dependency | `uv run pytest apps/api/tests/test_zpe_queue_targets.py` |
| Slurm VM control-plane config static validation | Tier-1 | Offline fixture-based deterministic checks | `scripts/validate-slurm-vm-offline.sh` |
| VM control-plane runtime health (`munge`, `slurmctld`, `slurmd`, `scontrol ping`, `sinfo`) | Tier-2 | Requires running Slurm services and host-level runtime state | `scripts/validate-slurm-vm-control-plane.sh --mode vm` |
| AiiDA -> Slurm minimal smoke (`verdi computer test core.slurm`) on VM | Tier-2 | Requires reachable Slurm controller and AiiDA scheduler integration | `scripts/aiida-slurm-smoke-vm.sh --artifact-dir investigations/artifacts/gra-89/manual` |
| Tier-2 one-command gate orchestration (control-plane + AiiDA smoke + summaries/bundle) | Tier-2 | Provides stable gate exit semantics and deterministic artifact packaging for release/cutover evidence | `scripts/aiida-tier2-gate-vm.sh --run-id <run-id>` |
| Remote rerun of VM smoke from lane machine | Tier-2 | Still real Slurm, but execution is remote via SSH | `scripts/aiida-slurm-smoke-vm-remote.sh --host <vm-host> --user <vm-user>` |
| Compatibility canary for AiiDA/Slurm/PostgreSQL/RabbitMQ/Python | Tier-2 (+ Tier-3 when out-of-range/major upgrades) | Runtime dependency safety must be validated before promotion/cutover | `docs/process/aiida-slurm-compatibility-matrix.md`<br>`scripts/aiida-vm-bootstrap.sh --env-file ops/aiida-vm/aiida-vm.env --sanity-check`<br>`scripts/aiida-slurm-smoke-vm.sh --artifact-dir investigations/artifacts/gra-132/<timestamp>/canary` |
| Full promotion-gate style remote VM flow (bootstrap + Slurm smoke) | Tier-3 | Broader operational confidence check; useful before cutover, not required per commit | `scripts/aiida-promotion-gate-vm-remote.sh --host <vm-host> --user <vm-user>` |
| Shared cluster validation of queue/account/qos policy in real environment | Tier-3 | Depends on external tenancy/policy and cluster availability; should be scheduled, not blocking | `scripts/aiida-slurm-smoke-vm.sh --profile <profile> --computer-label <shared-cluster-label>` |

## Practical Gate Recommendation

- PR gate default: Tier-1 only.
- Release/cutover gate: Tier-1 + Tier-2 (`scripts/aiida-tier2-gate-vm.sh` as the default Tier-2 entrypoint).
- Version/dependency upgrade gate: Tier-2 mandatory, Tier-3 mandatory for out-of-range or major upgrades (see `docs/process/aiida-slurm-compatibility-matrix.md`).
- Tier-3: run on demand (pre-cutover, incident repro, or policy drift checks).

## Notes on Existing E2E Tests

- `apps/api/tests/test_zpe_e2e_h2.py` and `apps/api/tests/test_zpe_e2e_molecules.py` require Quantum ESPRESSO binaries/pseudos, but they do not require Slurm.
- Keep them outside mandatory Slurm boundary checks unless intentionally executed on a Slurm-backed host for additional confidence.
