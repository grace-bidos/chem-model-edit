#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

OWNER="${RUNNER_OWNER:-grace-bidos}"
REPO="${RUNNER_REPO:-chem-model-edit}"
LABELS="${RUNNER_LABELS:-self-hosted,linux,x64,chem-trusted-pr}"

exec "${SCRIPT_DIR}/check_local_runner_health.sh" \
  --owner "${OWNER}" \
  --repo "${REPO}" \
  --labels "${LABELS}" \
  "$@"
