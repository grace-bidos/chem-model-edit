# GRA-83: Slurm Runtime Viability in Docker (Validation Scope)

Date: 2026-02-16  
Branch: `feature/gra-83-slurm-docker-viability`

## Objective

Assess whether Slurm can be used in Docker for validation in this project, with evidence from primary sources and lightweight local PoC checks.

## Bottom Line

Recommendation: **NO-GO as the primary validation runtime** (for production-like scheduler validation).  
Conditional allowance: **GO only for narrow smoke/integration checks** where we explicitly accept reduced fidelity.

Why:

- Slurm cgroup v2 behavior in Docker has explicit constraints that require elevated container privileges and host cgroup namespace alignment.
- Slurm daemons expect systemd/cgroup delegation behavior that is easy to diverge from in Docker test containers.
- Networking and hostname/port behavior must match Slurm controller/node assumptions; this is fragile in ad-hoc Docker topologies.
- MUNGE (or equivalent auth stack) setup and time-sync assumptions add extra failure modes in ephemeral containers.

## Concrete Constraints (Evidence-Based)

1. `systemd` and cgroup delegation
- Slurm cgroup v2 documentation requires proper systemd-managed daemon behavior (`Delegate=yes`) for cgroup controllers.
- Running daemons outside expected systemd delegation is documented as a special/limited path and not recommended for production equivalence.
- Source: Slurm cgroup v2 guide.

2. Docker privilege + cgroup namespace
- Slurm docs explicitly note Docker-specific requirements for cgroup v2 operation, including `--privileged` and `--cgroupns=host` in containerized daemon scenarios.
- Docker docs confirm `--privileged` broadens kernel capability/device access significantly.
- This materially reduces isolation and changes test risk posture.

3. Authentication (MUNGE / auth)
- Slurm quickstart documents MUNGE setup and key distribution requirements (plus time synchronization expectations).
- For transient Docker nodes/controllers, MUNGE key distribution and clock consistency can become a frequent source of false negatives.

4. Networking and name resolution
- `slurm.conf` documents explicit controller host/port settings (`SlurmctldHost`, `SlurmctldPort`, `SlurmdPort`) and node identity assumptions.
- Docker bridge/NAT/container DNS defaults can hide issues or create non-production behavior unless configured deliberately.

## Maintained Reference Stacks (Current Signals)

1. `giovtorres/slurm-docker-cluster`
- Purpose-fit for Docker-based Slurm experimentation and CI-friendly local validation.
- GitHub API snapshot (2026-02-16): `archived=false`, recent `pushed_at=2026-02-15T23:34:42Z`.
- Best use: local smoke and developer workflow checks.

2. `GoogleCloudPlatform/slurm-gcp`
- Actively maintained deployment stack for VM-based Slurm clusters on GCP (higher runtime fidelity than all-in-Docker daemons).
- GitHub API snapshot (2026-02-16): `archived=false`, `pushed_at=2026-02-09T14:52:37Z`.
- Best use: realistic validation baseline/fallback when scheduler correctness matters.

3. `NVIDIA/pyxis` (with Slurm container workflows)
- Maintained Slurm SPANK plugin for containerized job execution on Slurm clusters.
- GitHub API snapshot (2026-02-16): `archived=false`, `pushed_at=2026-02-12T06:06:49Z`.
- Best use: keep Slurm control plane on host/VM, containerize jobs instead of containerizing Slurm daemons.

## Lightweight Local PoC Checks (No Heavy Cluster Build)

Environment check results on this machine:

```bash
docker-found
29.1.3
CgroupDriver=cgroupfs CgroupVersion=1 SecurityOptions=["name=seccomp,profile=builtin"]
```

Container checks:

```bash
docker run --rm --cgroupns=private alpine:3.20 sh -c 'cat /proc/1/cgroup | head -n 3'
# returned cgroup entries successfully

docker run --rm --privileged --cgroupns=private alpine:3.20 sh -c 'test -w /sys/fs/cgroup && echo cgroup_writable'
# cgroup_writable
```

Interpretation:

- This host currently reports cgroup v1 via Docker info; cgroup-v2-specific Slurm behavior cannot be fully validated here.
- `--privileged` materially changes cgroup writability/capability context, supporting the documented risk that Docker-based validation may rely on elevated modes.

## Go/No-Go Decision and Fallback

Decision:

- **No-Go** for using Dockerized Slurm daemons as the authoritative validation runtime for scheduler behavior.
- **Go (limited)** for developer smoke tests and contract wiring checks only.

Fallback (recommended):

- Primary validation lane: run Slurm control-plane validation on VM/bare-metal-like environments (for example `slurm-gcp` style topology).
- Keep container usage at job-runtime layer (e.g., Pyxis/Enroot style) rather than daemon-runtime layer.
- If Docker must be used for fast feedback, gate it as non-blocking pre-checks and require a VM-based validation gate before merge/release.

## Sources (Primary, with Access Date)

Accessed: 2026-02-16

- Slurm cgroup v2 guide: https://slurm.schedmd.com/cgroup_v2.html
- Slurm quickstart (admin): https://slurm.schedmd.com/quickstart_admin.html
- Slurm configuration reference (`slurm.conf`): https://slurm.schedmd.com/slurm.conf.html
- Slurm containers guide: https://slurm.schedmd.com/containers.html
- Docker runtime privilege/capabilities: https://docs.docker.com/engine/containers/run/#runtime-privilege-and-linux-capabilities
- Docker `run` CLI reference (`--cgroupns`): https://docs.docker.com/reference/cli/docker/container/run/
- GitHub repository (slurm-docker-cluster): https://github.com/giovtorres/slurm-docker-cluster
- GitHub API (slurm-docker-cluster metadata): https://api.github.com/repos/giovtorres/slurm-docker-cluster
- GitHub repository (slurm-gcp): https://github.com/GoogleCloudPlatform/slurm-gcp
- GitHub API (slurm-gcp metadata): https://api.github.com/repos/GoogleCloudPlatform/slurm-gcp
- GitHub repository (pyxis): https://github.com/NVIDIA/pyxis
- GitHub API (pyxis metadata): https://api.github.com/repos/NVIDIA/pyxis
