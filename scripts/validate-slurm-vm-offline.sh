#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

"$ROOT_DIR/scripts/validate-slurm-vm-control-plane.sh" \
  --mode offline \
  --skip-runtime \
  --skip-key-perms
