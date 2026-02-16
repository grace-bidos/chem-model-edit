# Slurm VM Control-Plane Fixtures (GRA-87)

This directory contains deterministic fixtures used by
`scripts/validate-slurm-vm-control-plane.sh --mode offline`.

Files:

- `examples/slurm.conf`: minimal Slurm control-plane config for static validation.
- `examples/cgroup.conf`: minimal cgroup config for static validation.
- `examples/munge.key.example`: placeholder path for documentation only (never a real key).

Use cases:

- Validate parser/required-key checks in environments without a real VM.
- Keep validation deterministic when runtime services are intentionally unavailable.
