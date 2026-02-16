# WSL-to-VM Hardening Checklist for Slurm/AiiDA Validation

Use this checklist before running Slurm/AiiDA smoke and contract tests on a VM stack launched from WSL.

References (primary sources) are listed at the end with dates and links.

## 0) Preflight (must pass)

- `node`, `corepack`, `pnpm`, `uv` are available in the lane worktree.
- VM nodes are reachable over SSH from the operator environment.
- You can run `sudo` on controller and compute VMs.

## 1) WSL host baseline

### 1.1 Enable systemd in WSL distro

Set `/etc/wsl.conf` inside distro:

```ini
[boot]
systemd=true
```

Then restart WSL from Windows side:

```powershell
wsl --shutdown
```

Pass condition:

- `systemctl is-system-running` returns `running` or `degraded` (with understood non-critical units).

### 1.2 Standardize `.wslconfig` networking profile

Set in `%UserProfile%\.wslconfig`:

```ini
[wsl2]
networkingMode=mirrored
dnsTunneling=true
firewall=true
```

Restart WSL:

```powershell
wsl --shutdown
```

Pass condition:

- WSL distro has working DNS resolution and VM routes are reachable.

## 2) VM cgroup + systemd baseline

### 2.1 Require cgroup v2

Run on each VM:

```bash
stat -fc %T /sys/fs/cgroup
```

Pass condition:

- Output is `cgroup2fs`.

### 2.2 Configure Slurm cgroup policy

Set in `cgroup.conf` (controller policy distributed to nodes):

```ini
CgroupPlugin=cgroup/v2
ConstrainCores=yes
ConstrainRAMSpace=yes
ConstrainDevices=yes
```

Set in `slurm.conf` (controller policy distributed to nodes):

```ini
TaskPlugin=task/cgroup
```

Reload Slurm services after config distribution.

Pass condition:

- `scontrol show config | rg -i '^TaskPlugin\\s*=\\s*task/cgroup'`
- `scontrol show config | rg -i cgroup`
- `srun --mem=128M --cpus-per-task=1 --pty /bin/bash` succeeds and job-level limits are enforced.

### 2.3 Allow slurmd cgroup delegation

Create systemd override on compute nodes:

```bash
sudo systemctl edit slurmd
```

Add:

```ini
[Service]
Delegate=yes
```

Then:

```bash
sudo systemctl daemon-reload
sudo systemctl restart slurmd
```

Pass condition:

- `systemctl show slurmd -p Delegate` reports `Delegate=yes`.

## 3) Time synchronization gate (mandatory before MUNGE/Slurm)

### 3.1 Ensure chrony is active

```bash
sudo systemctl enable --now chronyd || sudo systemctl enable --now chrony
chronyc tracking
chronyc sources -v
```

Recommended `chrony.conf` line:

```ini
makestep 1.0 3
```

Pass condition:

- Tracking shows synchronized source and acceptable offset (team default target: < 100 ms).

### 3.2 Enforce MUNGE key security and clock sanity

```bash
ls -l /etc/munge/munge.key
```

Pass condition:

- Permissions are `-r--------` and owner/group is `munge munge` (or distro-default MUNGE service account if it differs).
- Clocks are synchronized on all controller/compute nodes before `munged` starts.

## 4) Networking and name-resolution gate

### 4.1 Validate Slurm control-plane naming

Check `slurm.conf` values (`SlurmctldHost`, `NodeName`, optional `NodeAddr`) and verify each VM can resolve all hosts:

```bash
getent hosts <controller-host>
getent hosts <compute-host>
```

Pass condition:

- Forward resolution succeeds on every participant VM.
- No participant resolves controller/node names to loopback addresses.

### 4.2 Restrict inbound firewall surface

Example with UFW on each VM:

```bash
sudo ufw default deny incoming
sudo ufw default allow outgoing
sudo ufw allow 22/tcp
sudo ufw allow 6817:6819/tcp
sudo ufw enable
sudo ufw status verbose
```

Pass condition:

- Only required management + Slurm ports are open.

## 5) SSH and MAC hardening gate

### 5.1 SSH daemon baseline

In `/etc/ssh/sshd_config` ensure:

```text
PasswordAuthentication no
PermitRootLogin no
PubkeyAuthentication yes
```

Then reload SSH:

```bash
sudo systemctl reload sshd || sudo systemctl reload ssh
```

Pass condition:

- Password login and root SSH login are blocked.

### 5.2 Keep AppArmor enabled

```bash
sudo aa-status
```

Pass condition:

- AppArmor is enabled; no global disablement for convenience.

## 6) Slurm + AiiDA validation gate

### 6.1 Slurm health

```bash
scontrol ping
sinfo
squeue
```

Pass condition:

- Controller responds and at least one node is schedulable.

### 6.2 AiiDA computer and scheduler path

```bash
verdi computer test <aiida-slurm-computer-label>
```

Then submit one minimal job through AiiDA to Slurm.

Pass condition:

- `verdi computer test` passes.
- AiiDA submission reaches `COMPLETED` with retrieved output.

## References (primary sources)

- Microsoft Learn, *Use systemd to manage Linux services with WSL* (2025-11-05): https://learn.microsoft.com/en-us/windows/wsl/systemd
- Microsoft Learn, *Advanced settings configuration in WSL* (2025-08-06): https://learn.microsoft.com/en-us/windows/wsl/wsl-config
- Linux kernel docs, *Control Group v2* (accessed 2026-02-16): https://www.kernel.org/doc/html/latest/admin-guide/cgroup-v2.html
- systemd docs, *systemd.resource-control(5), Delegate=* (accessed 2026-02-16): https://www.freedesktop.org/software/systemd/man/devel/systemd.resource-control.html
- SchedMD, *cgroup.conf* (page timestamp shows 2025-07-24): https://slurm.schedmd.com/cgroup.conf.html
- SchedMD, *cgroups guide* (accessed 2026-02-16): https://slurm.schedmd.com/cgroups.html
- SchedMD, *slurm.conf* (accessed 2026-02-16): https://slurm.schedmd.com/slurm.conf.html
- MUNGE Wiki, *Installation Guide* (2025-04-15): https://github.com/dun/munge/wiki/Installation-Guide
- chrony docs, *chrony.conf* (2024-09-02): https://chrony-project.org/doc/4.4/chrony.conf.html
- AiiDA docs, *How to run external codes on HPC resources* (accessed 2026-02-16): https://aiida.readthedocs.io/projects/aiida-core/en/stable/howto/run_codes.html
- AiiDA docs, *SLURM scheduler plugin API* (accessed 2026-02-16): https://aiida.readthedocs.io/projects/aiida-core/en/stable/reference/apidoc/aiida.schedulers.plugins.html#module-aiida.schedulers.plugins.slurm
- OpenSSH, *sshd_config(5)* (accessed 2026-02-16): https://man.openbsd.org/sshd_config
- Ubuntu docs, *UFW* (2025-10-17): https://help.ubuntu.com/community/UFW
- Ubuntu docs, *AppArmor* (2025-10-24): https://ubuntu.com/server/docs/security-apparmor
