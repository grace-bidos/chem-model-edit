# VM Baseline Provisioning from WSL (GRA-86)

This runbook provides a reproducible and low-risk baseline workflow for Hyper-V VM provisioning from WSL.

Scope:
- Create a baseline VM from WSL (safe defaults, idempotent bring-up)
- Capture a known-good snapshot
- Reset the VM back to that snapshot

## Files

- Entry point: `scripts/vm-baseline.sh`
- WSL wrapper: `ops/vm/vm-baseline.sh`
- Hyper-V implementation: `ops/vm/hyperv-baseline.ps1`

## Preflight

Run from repository root in WSL:

```bash
for tool in node corepack pnpm uv; do command -v "$tool" >/dev/null || { echo "missing: $tool"; exit 1; }; done
node -v
pnpm -v
uv --version
./scripts/vm-baseline.sh -Action preflight
```

The Hyper-V switch defaults to `Default Switch`. Override it with `-SwitchName` when needed.

## Minimal Bring-up (safe defaults)

```bash
./scripts/vm-baseline.sh -Action bringup
```

Default values:
- `VmName=chem-baseline`
- `SwitchName=Default Switch`
- `VhdPath=C:\HyperV\chem-baseline\chem-baseline.vhdx`
- `MemoryGb=8`
- `CpuCount=4`
- `DiskGb=120`

Notes:
- Bring-up is idempotent. If VM already exists, creation is skipped.
- VM is not auto-started unless `-StartVm` is passed.

## First-boot bring-up with installer ISO

```bash
./scripts/vm-baseline.sh -Action bringup \
  -IsoPath 'C:\ISO\ubuntu-24.04-live-server-amd64.iso' \
  -StartVm
```

## Baseline Snapshot Workflow

After initial OS setup and hardening in the VM, save the baseline checkpoint:

```bash
./scripts/vm-baseline.sh -Action snapshot -SnapshotName baseline
```

If you need to replace an existing baseline snapshot:

```bash
./scripts/vm-baseline.sh -Action snapshot -SnapshotName baseline -ReplaceSnapshot
```

## Reset Workflow

Reset VM state to the baseline snapshot:

```bash
./scripts/vm-baseline.sh -Action reset -SnapshotName baseline
```

Reset and start immediately:

```bash
./scripts/vm-baseline.sh -Action reset -SnapshotName baseline -StartVm
```

## Status

```bash
./scripts/vm-baseline.sh -Action status
```

Use status before and after reset to confirm snapshot-driven reproducibility.

## Safety Defaults

- No destructive action on bring-up for existing VMs.
- Snapshot replacement requires explicit `-ReplaceSnapshot`.
- Reset requires an existing snapshot name.
- VM auto-start is opt-in (`-StartVm`).
