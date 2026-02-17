#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Reconcile local ephemeral runner pool count.

Usage:
  scripts/runner/reconcile_ephemeral_pool.sh \
    --repo grace-bidos/chem-model-edit \
    --min 1 \
    --max 4 \
    --target 4

Options:
  --repo <owner/repo>     Required.
  --min <n>               Optional. Default: 1
  --baseline <n>          Deprecated alias for --min.
  --max <n>               Optional. Default: 4
  --target <n>            Optional. Default: max
  --lock-file <path>      Optional. Default: /run/lock/chem-model-edit-runner-pool.lock
  --dry-run               Print actions without mutating state.
  -h, --help              Show help.
EOF
}

repo=""
min_pool="1"
max_parallel="4"
target=""
lock_file="/run/lock/chem-model-edit-runner-pool.lock"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) repo="${2:-}"; shift 2 ;;
    --min) min_pool="${2:-}"; shift 2 ;;
    --baseline) min_pool="${2:-}"; shift 2 ;;
    --max) max_parallel="${2:-}"; shift 2 ;;
    --target) target="${2:-}"; shift 2 ;;
    --lock-file) lock_file="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$repo" ]]; then
  echo "--repo is required." >&2
  usage
  exit 1
fi

if ! command -v flock >/dev/null 2>&1; then
  echo "missing command: flock" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scale_script="${script_dir}/scale_local_runner_pool.sh"
if [[ ! -x "${scale_script}" ]]; then
  echo "missing executable: ${scale_script}" >&2
  exit 1
fi

if [[ -z "$target" ]]; then
  target="$max_parallel"
fi

cmd=(
  "${scale_script}"
  --repo "$repo"
  --min "$min_pool"
  --max "$max_parallel"
  --target "$target"
)
if [[ "$dry_run" -eq 1 ]]; then
  cmd+=(--dry-run)
fi

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] lock_file=${lock_file}"
  "${cmd[@]}"
  exit 0
fi

exec 9>"$lock_file"
if ! flock -n 9; then
  echo "Another reconcile process is running. lock=${lock_file}"
  exit 0
fi

"${cmd[@]}"
