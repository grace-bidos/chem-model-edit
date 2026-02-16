# Slurm VM Control-Plane Setup and Validation (GRA-87)

This runbook defines a practical setup/validation path for `munge`, `slurmctld`,
and `slurmd` on a single VM control-plane host.

## Goal

- Make VM prerequisites explicit before service startup.
- Validate config files deterministically, even when runtime execution is blocked.
- Provide a single helper command for operators.

## Artifacts

- Primary validator: `scripts/validate-slurm-vm-control-plane.sh`
- Offline deterministic wrapper: `scripts/validate-slurm-vm-offline.sh`
- Offline fixtures: `ops/slurm-vm/examples/*`

## Prerequisites (VM mode)

Required commands:

- `munge`
- `unmunge`
- `slurmctld`
- `slurmd`
- `scontrol`
- `sinfo`
- `systemctl` (recommended, for service state checks)

Required files (default paths):

- `/etc/slurm/slurm.conf`
- `/etc/slurm/cgroup.conf`
- `/etc/munge/munge.key`

Expected security baseline:

- `/etc/munge/munge.key` owner/group is `munge:munge`
- `/etc/munge/munge.key` mode is `0400`

## VM Setup Sequence

1. Install packages for Slurm and Munge on the VM.
2. Place a real Munge key at `/etc/munge/munge.key` with strict permissions.
3. Configure `/etc/slurm/slurm.conf` and `/etc/slurm/cgroup.conf`.
4. Enable and start services:
   - `sudo systemctl enable --now munge`
   - `sudo systemctl enable --now slurmctld`
   - `sudo systemctl enable --now slurmd`
5. Run validation:

```bash
scripts/validate-slurm-vm-control-plane.sh --mode vm
```

## Deterministic Validation Without a VM

Use this when running in CI, local laptops, or environments without Slurm services:

```bash
scripts/validate-slurm-vm-offline.sh
```

What this verifies:

- Required keys in `slurm.conf` and `cgroup.conf`
- `AuthType=auth/munge`
- At least one `NodeName` and one `PartitionName`

What this intentionally skips:

- Service state checks
- Runtime checks (`munge -n | unmunge`, `SLURM_CONF=<path> slurmctld -V`, `SLURM_CONF=<path> slurmd -C`, `SLURM_CONF=<path> scontrol ping`)

## Exit Codes

- `0`: all checks passed
- `1`: deterministic config failures
- `2`: runtime/prerequisite blockers (static checks still executed)

## Typical Blockers

Examples of actionable blockers reported by the validator:

- missing required command (`slurmctld`, `munge`, etc.)
- missing config/key path (`/etc/slurm/slurm.conf`, `/etc/munge/munge.key`)
- service not active (`munge`, `slurmctld`, `slurmd`)
- controller unreachable (`scontrol ping` failure)

## Observed Blockers in This Development Environment

Execution date: **2026-02-16**

Command:

```bash
scripts/validate-slurm-vm-control-plane.sh --mode vm
```

Observed result:

- Exit code: `2`
- Blockers:
  - `missing slurm.conf: /etc/slurm/slurm.conf`
  - `missing cgroup.conf: /etc/slurm/cgroup.conf`
  - `missing munge key: /etc/munge/munge.key`

Interpretation: this workspace host is not configured as a real Slurm VM control-plane,
so runtime checks are blocked by missing host-level configuration files.

## Suggested Operator Loop

1. Run `--mode vm`.
2. Fix all blockers and deterministic failures.
3. Re-run until exit code is `0`.
