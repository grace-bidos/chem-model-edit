# GRA-89 AiiDA->Slurm VM Smoke Investigation (2026-02-16)

## Scope

- Add and validate a VM-oriented smoke execution path for AiiDA -> Slurm.
- Ensure blocked-path evidence collection and fallback diagnostics are executable.

## Preflight

Executed in `/home/grace/projects/chem-model-edit/.worktrees/gra-89`:

- `for tool in node corepack pnpm uv; do command -v "$tool" >/dev/null || exit 1; done` -> pass
- `node -v` -> `v22.16.0`
- `pnpm -v` -> `10.27.0`
- `uv --version` -> `uv 0.6.14`

## Added Assets

- Script: `scripts/aiida-slurm-smoke-vm.sh`
- Runbook: `docs/process/aiida-slurm-smoke-vm.md`

## Validation Performed

Syntax/help checks:

- `bash -n scripts/aiida-slurm-smoke-vm.sh` -> pass
- `scripts/aiida-slurm-smoke-vm.sh --help` -> pass

Smoke execution in this environment:

```bash
scripts/aiida-slurm-smoke-vm.sh
```

Result:

- Exit code: `2` (blocked path expected in this non-VM environment)

Artifacts:

- `investigations/artifacts/gra-89/20260216-093715-vm-slurm/`

Key evidence from logs:

- `diag-which-verdi.log`: `verdi` not found on PATH.
- `diag-slurm-scontrol-ping.log`: Slurm configuration source unavailable.
- `diag-slurm-sinfo.log`: `Could not establish a configuration source`.
- `diag-slurm-squeue-user.log`: same Slurm configuration blocker.
- `fallback-next-steps.txt`: actionable VM-side recovery checklist.

## Blocker Interpretation

This host does not satisfy VM smoke prerequisites:

- Missing AiiDA CLI (`verdi`).
- Slurm binaries exist but controller/config source is unresolved.

The script correctly classified this as `blocked` and emitted fallback diagnostics for VM execution.

## VM Re-run Command

Run on the target VM with AiiDA + Slurm configured:

```bash
scripts/aiida-slurm-smoke-vm.sh --profile <aiida_profile_name>
```

Expected result on healthy VM:

- Exit code `0`
- `aiida-computer-test-slurm.log` shows scheduler smoke success.
