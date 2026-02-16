# Remote VM Validation Driver for Slurm+AiiDA Smoke Reruns (GRA-124)

This runbook defines the remote rerun driver that executes VM smoke checks over SSH and collects evidence artifacts locally with deterministic naming.

## Goal

- Re-run `scripts/aiida-slurm-smoke-vm.sh` on a remote VM from the local lane environment.
- Keep remote execution simple (`ssh` only) while preserving blocker semantics.
- Store local evidence in a deterministic path for each rerun.

## Script

- `scripts/aiida-slurm-smoke-vm-remote.sh`

## Required Environment Variables

- `REMOTE_VM_HOST`: target VM host (`hostname`, IP, or `user@host`)
- `REMOTE_VM_REPO_DIR`: repository root path on the remote VM where smoke script exists

## Optional Environment Variables

- `REMOTE_VM_USER`: SSH user if not included in `REMOTE_VM_HOST`
- `REMOTE_VM_PORT`: SSH port (default: `22`)
- `REMOTE_VM_SMOKE_SCRIPT`: remote smoke script path (default: `scripts/aiida-slurm-smoke-vm.sh`)
- `REMOTE_VM_SSH_KEY`: SSH key path
- `REMOTE_SSH_BIN`: SSH binary path (default: `ssh`)
- `LOCAL_ARTIFACT_BASE`: local evidence root (default: `investigations/artifacts/gra-124`)
- `REMOTE_ARTIFACT_BASE`: remote evidence root (default: `/tmp/gra-124-remote-vm-validation`)
- `GRA124_RUN_ID`: deterministic run id override (`YYYYMMDDTHHMMSSZ` format recommended)
- `AIIDA_PROFILE`: forwarded to remote smoke script (`--profile`)
- `AIIDA_COMPUTER_LABEL`: forwarded to remote smoke script (`--computer-label`)
- `AIIDA_WORKDIR_TEMPLATE`: forwarded to remote smoke script (`--workdir`)

## Windows OpenSSH Path Compatibility

The driver accepts Windows-style paths for `REMOTE_SSH_BIN` and `REMOTE_VM_SSH_KEY` and normalizes them for WSL execution.

Examples:

- `C:\Windows\System32\OpenSSH\ssh.exe`
- `C:\Users\grace\.ssh\id_ed25519`

Equivalent explicit WSL paths also work:

- `/mnt/c/Windows/System32/OpenSSH/ssh.exe`
- `/mnt/c/Users/grace/.ssh/id_ed25519`

## Usage

Minimal invocation:

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
scripts/aiida-slurm-smoke-vm-remote.sh
```

Explicit rerun id and profile:

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
GRA124_RUN_ID=20260216T120000Z \
AIIDA_PROFILE=gra89-vm \
scripts/aiida-slurm-smoke-vm-remote.sh
```

Windows OpenSSH path through WSL:

```bash
REMOTE_VM_HOST=Administrator@10.0.0.40 \
REMOTE_VM_REPO_DIR=/home/admin/chem-model-edit \
REMOTE_SSH_BIN='C:\Windows\System32\OpenSSH\ssh.exe' \
REMOTE_VM_SSH_KEY='C:\Users\grace\.ssh\id_ed25519' \
scripts/aiida-slurm-smoke-vm-remote.sh --ssh-opt StrictHostKeyChecking=accept-new
```

Dry-run planning only:

```bash
REMOTE_VM_HOST=operator@vm-controller \
REMOTE_VM_REPO_DIR=/home/operator/chem-model-edit \
scripts/aiida-slurm-smoke-vm-remote.sh --dry-run
```

## Exit Codes

- `0`: remote smoke passed and local evidence collection succeeded
- `1`: local invocation/config error or hard remote script failure
- `2`: blocked (SSH/precondition blocker, remote smoke blocker, or evidence collection blocker)

## Evidence Artifacts

Local run directory format:

- `<LOCAL_ARTIFACT_BASE>/<target-slug>/<run-id>-<target-slug>-slurm-aiida-smoke/`

Typical evidence files:

- `run-metadata.txt`
- `remote-preflight.log`
- `remote-smoke.log`
- `remote-artifacts.tar.gz`
- `remote-artifacts/` (expanded archive)
- `fallback-next-steps.txt` (on blocked/hard-failure path)

The naming scheme is deterministic for a given `GRA124_RUN_ID` and target.

