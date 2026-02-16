# AiiDA->Slurm VM Smoke Execution (GRA-89)

This runbook defines a VM-focused smoke path for validating the minimal AiiDA -> Slurm submission chain.

## Goal

- Confirm the VM has usable Slurm control-plane connectivity.
- Confirm AiiDA profile availability on the VM.
- Confirm AiiDA can perform scheduler-level smoke against `core.slurm` via `verdi computer test`.
- If blocked, capture evidence and provide immediate fallback diagnostics.

## Preconditions

Required commands on the VM:

- `verdi`
- `sinfo`
- `sbatch`
- `squeue`
- `scontrol`

Operational assumptions:

- Slurm controller (`slurmctld`) is reachable from the VM.
- At least one AiiDA profile is configured on the VM.
- VM user has permissions for the target Slurm partition/account policy.

Known host-side caveats before VM smoke:

- Avoid VM address plans in `172.17.0.0/16` to prevent route conflicts from WSL/container bridge defaults.
- If Windows OpenSSH is in the operator path, use canonical key path and strict ACLs (`C:\\ProgramData\\ssh\\administrators_authorized_keys`).
- If Secure Boot DB hash mismatch blocks boot/runtime, align image signing with firmware template or apply a documented temporary disable/re-enable workaround.

## Script

- `scripts/aiida-slurm-smoke-vm.sh`

Default behavior:

1. Captures environment metadata (`hostname`, `uname`, `whoami`).
2. Verifies runtime command availability.
3. Resolves AiiDA profile (explicit `--profile` or current default profile).
4. Runs Slurm control checks (`scontrol ping`, `sinfo`).
5. Creates/reuses AiiDA `core.slurm` computer (`gra89-vm-slurm`).
6. Runs `verdi -p <profile> computer test <computer>`.
7. On blocked path, collects fallback diagnostics and writes next steps.

## Commands

Default run (uses current default AiiDA profile):

```bash
scripts/aiida-slurm-smoke-vm.sh
```

Run against explicit profile:

```bash
scripts/aiida-slurm-smoke-vm.sh --profile <aiida_profile_name>
```

Override computer label/workdir:

```bash
scripts/aiida-slurm-smoke-vm.sh \
  --profile <aiida_profile_name> \
  --computer-label gra89-vm-slurm-custom \
  --workdir /tmp/aiida-gra89-custom-{username}
```

Override artifact directory:

```bash
scripts/aiida-slurm-smoke-vm.sh \
  --artifact-dir investigations/artifacts/gra-89-vm
```

## Exit Codes

- `0`: smoke passed
- `1`: hard failure (invalid invocation/unexpected script failure)
- `2`: blocked (missing runtime/profile/Slurm reachability, or AiiDA Slurm computer test failure)

## Artifacts

Logs are written under:

- `investigations/artifacts/gra-89/<timestamp>-vm-slurm/`

Typical files:

- `slurm-scontrol-ping.log`
- `slurm-sinfo.log`
- `aiida-computer-test-slurm.log`
- `diag-verdi-profile-list.log` (when blocked)
- `fallback-next-steps.txt` (when blocked)

Minimum accepted evidence set for a successful run (`exit 0`):

- `slurm-scontrol-ping.log` shows controller `UP`.
- `slurm-sinfo.log` shows at least one schedulable node.
- `aiida-computer-test-slurm.log` ends with all checks succeeded.
- `aiida-computer-setup-slurm.log` when this run created the computer, or a short operator note indicating expected reuse for this label.
- Runtime identity/context logs exist (`env-hostname.log`, `env-uname.log`, `env-whoami.log`).

This aligns with the GRA-94 VM validation runbook evidence baseline.

Optional but recommended:

- `env-verdi-version.log`
- `env-slurm-version.log`
- a short operator note recording Secure Boot status and any temporary exception used.

## Blocked Path Interpretation

Common blockers include:

- Missing `verdi` on VM path.
- No resolvable AiiDA profile.
- Slurm controller not reachable (`scontrol ping`/`sinfo` failure).
- AiiDA `computer test` scheduler failure due to account/partition/policy mismatch.

When blocked, inspect logs in the run artifact directory first, then execute the suggested checks in `fallback-next-steps.txt`.
