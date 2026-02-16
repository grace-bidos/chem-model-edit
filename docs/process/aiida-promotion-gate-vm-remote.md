# Remote VM Promotion Gate Rerun Driver (GRA-126)

This runbook defines the remote rerun driver for the PostgreSQL + RabbitMQ promotion gate VM checks.
The driver executes gate 3/4 from `docs/process/aiida-promotion-gate-vm.md` over SSH and collects evidence locally.

## Goal

- Re-run gate 3 (VM bootstrap verification) and gate 4 (AiiDA -> Slurm VM smoke) on a remote VM.
- Keep execution transport simple (`ssh` only).
- Persist deterministic local evidence for the promotion decision record.

## Script

- `scripts/aiida-promotion-gate-vm-remote.sh`

## Required Environment Variables

- `REMOTE_VM_HOST`: target VM host (`hostname`, IP, or `user@host`)
- `REMOTE_VM_REPO_DIR`: repository root path on the remote VM

## Optional Environment Variables

- `REMOTE_VM_USER`: SSH user when not included in `REMOTE_VM_HOST`
- `REMOTE_VM_PORT`: SSH port (default: `22`)
- `REMOTE_VM_BOOTSTRAP_SCRIPT`: remote bootstrap script (default: `scripts/aiida-vm-bootstrap.sh`)
- `REMOTE_VM_SMOKE_SCRIPT`: remote smoke script (default: `scripts/aiida-slurm-smoke-vm.sh`)
- `REMOTE_AIIDA_ENV_FILE`: remote env file for bootstrap (default: `ops/aiida-vm/aiida-vm.env`)
- `REMOTE_VM_SSH_KEY`: SSH private key path
- `REMOTE_SSH_BIN`: SSH binary path (default: `ssh`)
- `LOCAL_ARTIFACT_BASE`: local evidence root (default: `investigations/artifacts/gra-126`)
- `REMOTE_ARTIFACT_BASE`: remote evidence root (default: `/tmp/gra-126-promotion-run-vm`)
- `GRA126_RUN_ID`: deterministic run id override (`YYYYMMDDTHHMMSSZ` recommended)
- `AIIDA_PROFILE`: forwarded to slurm smoke script (`--profile`)
- `AIIDA_COMPUTER_LABEL`: forwarded to slurm smoke script (`--computer-label`)
- `AIIDA_WORKDIR_TEMPLATE`: forwarded to slurm smoke script (`--workdir`)

## Usage

Remote sanity-check path (no `--apply`):

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
scripts/aiida-promotion-gate-vm-remote.sh
```

First-time bootstrap path (explicit apply/init-profile):

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
REMOTE_AIIDA_ENV_FILE=ops/aiida-vm/aiida-vm.env \
scripts/aiida-promotion-gate-vm-remote.sh \
  --bootstrap-apply \
  --bootstrap-init-profile
```

Deterministic rerun id and explicit profile:

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
GRA126_RUN_ID=20260216T120000Z \
AIIDA_PROFILE=chem-model-vm \
scripts/aiida-promotion-gate-vm-remote.sh
```

Dry-run planning only:

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
scripts/aiida-promotion-gate-vm-remote.sh --dry-run
```

## Exit Codes and Blocker Interpretation

- `0`: gate 3/4 remote run passed and artifact collection succeeded (`PROMOTE`-eligible from VM runtime perspective)
- `1`: hard failure in invocation or unexpected remote failure (fix script/operator error before rerun)
- `2`: blocked path (`HOLD`) due to one of:
  - SSH/preflight failure
  - gate 3 bootstrap/verdi checks failed
  - gate 4 slurm smoke reported blocked state
  - remote evidence collection failed

When blocked, use:

- `fallback-next-steps.txt` for immediate rerun commands
- `promotion-gate-vm-summary.md` for final exit interpretation

## Evidence Artifacts

Local run directory format:

- `<LOCAL_ARTIFACT_BASE>/<target-slug>/<run-id>-<target-slug>-promotion-gate-vm/`

Typical files:

- `run-metadata.txt`
- `remote-preflight.log`
- `remote-promotion-run.log`
- `promotion-gate-vm-summary.md`
- `remote-artifacts.tar.gz`
- `remote-artifacts/<run-name>/03-vm-bootstrap.log`
- `remote-artifacts/<run-name>/03-verdi-profile-list.log`
- `remote-artifacts/<run-name>/03-verdi-status.log`
- `remote-artifacts/<run-name>/04-slurm-smoke-driver.log`
- `remote-artifacts/<run-name>/04-slurm-smoke-summary.md`
- `fallback-next-steps.txt` (blocked/hard-failure path)
