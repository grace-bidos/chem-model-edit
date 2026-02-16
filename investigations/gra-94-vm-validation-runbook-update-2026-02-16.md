# GRA-94 VM Validation Runbook Update Findings (2026-02-16)

## Scope

Documented VM validation findings into runbooks for lane-B documentation-only delivery:

- WSL route caveat for `172.17.x.x` ranges
- Windows OpenSSH path/ACL requirements
- Secure Boot DB hash mismatch issue and workaround policy
- minimum accepted evidence set for successful AiiDA -> Slurm VM smoke

## Updated runbooks

- `docs/process/wsl-to-vm-hardening-checklist.md`
- `docs/process/aiida-slurm-smoke-vm.md`

## Findings captured

### 1) WSL routing caveat to `172.17.x.x`

Using `172.17.0.0/16` for VM addressing can collide with container bridge defaults and misroute traffic from the operator environment.

Operational policy:

- avoid `172.17.0.0/16` for Slurm/AiiDA VM subnets
- use a dedicated non-conflicting subnet, or explicitly pin routes as a temporary control

### 2) Windows OpenSSH path/ACL requirements

When Windows OpenSSH is part of the operator path, admin key auth should use:

- `C:\\ProgramData\\ssh\\administrators_authorized_keys`

Permissions should remain tightly scoped to admin/system principals to avoid OpenSSH permission rejection and accidental key exposure.

### 3) Secure Boot DB hash issue

On some Generation 2 VM/image combinations, signed binary verification can fail due to firmware DB hash/certificate mismatch.

Operational policy:

- preferred: align image/kernel/shim signing chain with firmware secure boot template
- temporary exception: disable secure boot only for bootstrap/validation, then re-enable after alignment

### 4) Minimum accepted evidence set from successful run

A successful smoke (`exit 0`) should include, at minimum:

- `slurm-scontrol-ping.log`
- `slurm-sinfo.log`
- `aiida-computer-test-slurm.log`
- `aiida-computer-setup-slurm.log` (or explicit reuse evidence)
- `env-hostname.log`, `env-uname.log`, `env-whoami.log`

Reference successful artifact example in repository root worktree:

- `investigations/artifacts/gra-89-vm/20260216-052824-vm-slurm/`
