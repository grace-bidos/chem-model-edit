# GRA-90 WSL-to-VM Hardening Research for Slurm/AiiDA Validation (2026-02-16)

## Scope

Research objective: produce a primary-source-backed hardening baseline for a Windows/WSL host running Linux VMs used for Slurm + AiiDA validation.

Focus areas:

- cgroup and systemd behavior for Slurm job isolation
- time synchronization requirements (critical for MUNGE and scheduler stability)
- networking model and host resolution
- security controls for SSH, MUNGE key material, and host firewalling

## Primary-Source Findings (dated)

| Area | Primary source | Date on source | Key finding | Concrete recommendation |
| --- | --- | --- | --- | --- |
| WSL init/systemd | Microsoft Learn: [Use systemd to manage Linux services with WSL](https://learn.microsoft.com/en-us/windows/wsl/systemd) | 2025-11-05 | WSL supports systemd via `/etc/wsl.conf` (`[boot] systemd=true`). | Require systemd-enabled distro before any Slurm service bring-up. |
| WSL networking | Microsoft Learn: [WSL config (`.wslconfig`)](https://learn.microsoft.com/en-us/windows/wsl/wsl-config) | 2025-08-06 | `networkingMode=mirrored` is available, with `dnsTunneling=true` and `firewall=true` defaults documented. | Standardize mirrored networking in hardened profile; keep DNS tunneling and firewall integration enabled unless an explicit exception is approved. |
| Linux cgroup v2 | Linux kernel docs: [Control Group v2](https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html) | (no explicit page date shown, accessed 2026-02-16) | Unified hierarchy semantics and delegation constraints define safe resource-control boundaries. | Validate VM nodes run unified cgroup v2 hierarchy before enabling Slurm cgroup enforcement. |
| systemd + cgroup delegation | `systemd.resource-control(5)`: [Delegate=](https://www.freedesktop.org/software/systemd/man/devel/systemd.resource-control.html) | (devel docs, accessed 2026-02-16) | `Delegate=` is required when a service manages cgroups below its unit. | Add `Delegate=yes` override to `slurmd.service` on cgroup-v2 nodes. |
| Slurm cgroup plugin | SchedMD: [cgroup.conf](https://slurm.schedmd.com/cgroup.conf.html) and [cgroups guide](https://slurm.schedmd.com/cgroups.html) | cgroup.conf page shows 2025-07-24 | Slurm recommends cgroup/v2 path and provides constraints knobs (`ConstrainCores`, `ConstrainRAMSpace`, `ConstrainDevices`, etc.). | Pin cgroup plugin policy explicitly and enable CPU/memory/device constraints for validation jobs. |
| Slurm host identity | SchedMD: [slurm.conf](https://slurm.schedmd.com/slurm.conf.html) | (page current as accessed 2026-02-16) | `SlurmctldHost` and node naming are central to controller/node communication. | Use stable FQDN/IP mapping for controller and compute VMs; verify forward+reverse resolution before smoke tests. |
| MUNGE clock dependency + key perms | MUNGE wiki: [Installation Guide](https://github.com/dun/munge/wiki/Installation-Guide) | 2025-04-15 | Credentials are time-bound; docs warn clocks must be synchronized. Key should be owner-only readable. | Enforce NTP sync gate and keep `/etc/munge/munge.key` at `0400` owned by `munged:munged`. |
| Time sync mechanism | chrony docs: [chrony.conf](https://chrony-project.org/doc/4.4/chrony.conf.html) | 2024-09-02 | `makestep` allows fast correction of large startup offsets. | Prefer `chronyd` in VMs, set bounded `makestep`, and fail validation if drift exceeds threshold. |
| AiiDA + Slurm integration | AiiDA docs: [How to run external codes on HPC resources](https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/run_codes.html) and [SLURM scheduler plugin](https://aiida.readthedocs.io/projects/aiida-core/en/stable/reference/apidoc/aiida.schedulers.plugins.html#module-aiida.schedulers.plugins.slurm) | (stable docs, accessed 2026-02-16) | AiiDA requires correctly configured computer transport + scheduler metadata for remote execution. | Gate acceptance on `verdi computer test` success plus a real Slurm submission/parse cycle. |
| SSH hardening controls | OpenSSH man page: [sshd_config(5)](https://man.openbsd.org/sshd_config) | (manpage current as accessed 2026-02-16) | Security-relevant controls include disabling password auth and root login. | Use key-only auth and disable root/password login on Slurm/AiiDA VMs. |
| Host firewalling | Ubuntu docs: [UncomplicatedFirewall](https://help.ubuntu.com/community/UFW) | 2025-10-17 | UFW enables default-deny host policy with explicit allow rules. | Default deny incoming; explicitly allow SSH and Slurm ports only. |
| Mandatory access control | Ubuntu docs: [AppArmor](https://ubuntu.com/server/docs/security-apparmor) | 2025-10-24 | AppArmor should remain enabled with enforce mode where possible. | Keep AppArmor enabled; do not disable globally for bootstrap convenience. |

## Recommended Hardening Profile (WSL host -> Linux VM)

### 1) cgroup/systemd baseline

- Require systemd-enabled WSL distro (`/etc/wsl.conf` with `[boot] systemd=true`).
- Require VM kernel in unified cgroup v2 mode.
- Set Slurm cgroup policy explicitly:
  - `CgroupPlugin=cgroup/v2` (or controlled `autodetect` with cgroup-v2 validation)
  - `ConstrainCores=yes`
  - `ConstrainRAMSpace=yes`
  - `ConstrainDevices=yes`
- Add `slurmd.service` override with `Delegate=yes`.

### 2) time sync gates

- Run `chronyd` on every Slurm/AiiDA VM.
- Use `makestep` in `chrony.conf` for startup correction.
- Fail fast if time is unsynchronized before starting `munged`, `slurmctld`, `slurmd`, or AiiDA daemons.

### 3) networking gates

- Prefer WSL mirrored networking profile for predictable host/VM routing.
- Keep DNS tunneling enabled unless there is a documented resolver conflict.
- Validate name resolution for `SlurmctldHost` and compute nodes from every participant VM.

### 4) security gates

- SSH: `PasswordAuthentication no`, `PermitRootLogin no`, key-based auth only.
- MUNGE key: `0400` and owned by `munged:munged`.
- Host firewall: default deny incoming, then allow minimal required ports.
- Keep AppArmor enabled and investigate policy exceptions instead of disabling enforcement.

## Acceptance criteria for GRA-90 output

A hardened lane is ready for downstream VM smoke validation when all checks pass:

- systemd enabled in WSL guest distro and Slurm services managed by systemd.
- cgroup v2 active on VM nodes and Slurm cgroup constraints enabled.
- NTP synchronized state confirmed before MUNGE/Slurm startup.
- Controller/node hostname resolution is stable from all VMs.
- SSH + firewall + MUNGE key permissions satisfy minimum hardening policy.
- AiiDA `verdi computer test` passes and one Slurm-backed AiiDA test job completes.
