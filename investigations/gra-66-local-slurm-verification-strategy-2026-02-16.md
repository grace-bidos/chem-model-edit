# GRA-66A: Reproducible Local/Dev Slurm Verification Strategy (2026-02-16)

## Objective

Define a reproducible verification strategy for local/dev Slurm validation that is explicit about fidelity limits, setup cost, and CI suitability for this repository.

## Executive Recommendation

Use a **tiered strategy**:

1. **Primary local verification lane: VM with real Slurm daemons (recommended baseline).**
2. **Secondary fast lane: mock/contract tests** for API and scheduler-adapter behavior.
3. **Optional developer smoke lane: containerized Slurm** for non-blocking checks only.
4. **Promotion/escape hatch: remote shared cluster** for high-fidelity confirmation when local environment is constrained.

This ordering reflects hard constraints in Slurm cgroup v2 container operation and AiiDA quick-install limitations.

## Primary-Source Anchors

### Slurm constraints that affect reproducibility

- Slurm quickstart requires synchronized clocks/users/groups, consistent auth material (`munge.key`), and daemon prerequisites across nodes ([Slurm Quick Start Admin](https://slurm.schedmd.com/quickstart_admin.html)).
- A minimal operational cluster still requires at least controller + compute role and correct service/user provisioning ([Slurm Quick Start Admin](https://slurm.schedmd.com/quickstart_admin.html)).
- `cgroup.conf` marks several options as development/testing oriented and warns about unsafe combinations on systemd systems (`IgnoreSystemd`, `CgroupMountpoint`) ([cgroup.conf](https://slurm.schedmd.com/cgroup.conf.html)).
- Slurm cgroup v2 docs state container namespace/mount requirements and Docker caveats: write access to cgroups requires `--privileged`, and host namespace mode needs `--cgroupns=host` with `--cgroup-parent` ([cgroup_v2](https://slurm.schedmd.com/cgroup_v2.html)).

### AiiDA constraints that affect verification scope

- AiiDA quick install (`verdi presto`) defaults to SQLite and no RabbitMQ unless services are available; docs explicitly call out functionality/performance limits ([Quick install guide](https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_quick.html)).
- Without RabbitMQ, quick setup cannot run daemon-managed flow (`submit`, daemon control, play/pause/kill process controls) ([Quick install guide](https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_quick.html)).
- Complete install positions RabbitMQ as required for daemon execution and `core.psql_dos` as the production-performance storage path ([Complete install guide](https://aiida.readthedocs.io/projects/aiida-core/en/stable/installation/guide_complete.html)).
- Computer onboarding for scheduler-backed execution is a two-step setup/configure flow (`verdi computer setup`, `verdi computer configure ...`) and should be validated with `verdi computer test` ([How to run external codes](https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/run_codes.html)).

## Recommended Strategy by Verification Layer

### Layer 1 (default gate): VM real Slurm + AiiDA service profile

**Why this is the baseline**

- Closest match to Slurm daemon assumptions (service user, systemd/cgroup lifecycle, munge/key distribution semantics) from quickstart and cgroup docs.
- Avoids Docker-specific cgroup v2 privilege/namespace caveats that reduce equivalence.
- Supports realistic AiiDA flow with PostgreSQL + RabbitMQ + daemon.

**Setup steps (reproducible local/dev)**

1. Provision a Linux VM (systemd enabled, cgroup v2 capable).
2. Install and configure Slurm per quickstart: synced clocks, users/groups, shared `munge.key`, `slurm` user/service directories, controller + compute roles.
3. Configure `cgroup.conf` conservatively:
   - `CgroupPlugin=autodetect` or explicit `cgroup/v2`.
   - Avoid `IgnoreSystemd*` except temporary diagnostics.
4. Start and validate Slurm daemons (`slurmctld`, `slurmd`) and basic scheduling commands.
5. Install AiiDA with complete profile path:
   - `core.psql_dos` storage.
   - RabbitMQ enabled.
6. Register compute target in AiiDA:
   - `verdi computer setup`
   - `verdi computer configure core.local|core.ssh <label>`
   - `verdi computer test <label>`
7. Run smoke workflow that includes submission path (daemon-backed), not only direct/local execution.

**Explicit limitations**

- Higher setup/runtime cost than mocks or containers.
- Requires VM management discipline (image pinning, repeatable bootstrap scripts).

### Layer 2 (always-on fast checks): mock/contract tests

**Why include**

- Deterministic and cheap; ideal for PR-level confidence on adapter contracts and payload mapping.
- Does not depend on scheduler service health.

**Setup steps**

1. Define stable scheduler-facing contracts (submit, poll, cancel, parse status).
2. Add fixture-based responses for success, pending, failure, timeout, and malformed outputs.
3. Run in CI as required checks.

**Explicit limitations**

- Cannot validate real Slurm daemon behavior, cgroup enforcement, or AiiDA daemon/broker integration.

### Layer 3 (optional non-blocking): containerized Slurm smoke

**Why include**

- Fast local feedback when VM is unavailable.

**Setup steps**

1. Use documented cgroup v2 container requirements:
   - ensure valid cgroup/mount/process namespaces,
   - for Docker use either private cgroup namespace or host mode with parent cgroup,
   - account for `--privileged` requirement for writable cgroup operations.
2. Run reduced-scope smoke only (controller reachability/basic command path).

**Explicit limitations**

- Docker privilege/namespace requirements reduce production equivalence and can mask or create behavior not seen on VM/bare metal.
- Treat as informational/non-blocking.

### Layer 4 (promotion fallback): remote shared cluster

**Why include**

- High-fidelity confirmation when local host constraints block reliable Slurm bring-up.

**Setup steps**

1. Configure AiiDA remote computer via SSH transport.
2. Validate remote scheduler connectivity with `verdi computer test`.
3. Execute bounded smoke runs with explicit queue/resource limits.

**Explicit limitations**

- Lower determinism and higher operational coupling (queue contention, policy differences, network/access dependencies).

## Reproducibility Controls (required regardless of lane)

- Pin versions for Slurm, AiiDA, PostgreSQL, RabbitMQ in runbooks.
- Keep bootstrap scripts idempotent and log artifacts under `investigations/artifacts/`.
- Separate blocking vs non-blocking gates:
  - Blocking: mock/contract + VM real-Slurm smoke.
  - Non-blocking: containerized Slurm smoke.
- Track environment fingerprint in every run (kernel, cgroup mode, Slurm/AiiDA versions).

## Decision Matrix

| Option | Fidelity to production scheduler behavior | Reproducibility | Setup/maintenance cost | CI speed | Recommended role |
|---|---|---|---|---|---|
| VM real Slurm | High | High (if image/bootstrap pinned) | Medium-High | Medium | **Primary local/dev gate** |
| Containerized Slurm | Low-Medium (cgroup/namespace caveats) | Medium | Medium | Fast | Optional, non-blocking smoke |
| Mock/contract tests | Low (runtime) / High (API contract) | Very High | Low | Very Fast | Always-on required checks |
| Remote shared cluster | High | Medium-Low (shared variability) | Medium | Slow-Medium | Promotion fallback / confirmation |

## Final Recommendation

Adopt **VM real Slurm + daemon-capable AiiDA profile** as the canonical local/dev verification path, keep **mock/contract tests** as mandatory fast checks, and treat **containerized Slurm** as advisory only. Use **remote shared cluster** as fallback when local VM constraints prevent timely verification.
