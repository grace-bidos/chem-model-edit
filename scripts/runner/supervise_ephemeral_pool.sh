#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Supervise local ephemeral runner pool in a continuous loop.

Usage:
  scripts/runner/supervise_ephemeral_pool.sh \
    --repo grace-bidos/chem-model-edit \
    --min 1 \
    --max 4 \
    --target 4 \
    --interval 15

Options:
  --repo <owner/repo>       Required.
  --min <n>                 Optional. Default: 1
  --baseline <n>            Deprecated alias for --min.
  --max <n>                 Optional. Default: 4
  --target <n>              Optional. Default: max
  --interval <seconds>      Optional. Default: 15
  --lock-file <path>        Optional. Default: /run/lock/chem-model-edit-runner-pool.lock
  --once                    Run one reconcile iteration and exit.
  --dry-run                 Print actions only.
  -h, --help                Show help.
EOF
}

repo=""
min_pool="1"
max_parallel="4"
target=""
interval_seconds="15"
lock_file="/run/lock/chem-model-edit-runner-pool.lock"
once=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) repo="${2:-}"; shift 2 ;;
    --min) min_pool="${2:-}"; shift 2 ;;
    --baseline) min_pool="${2:-}"; shift 2 ;;
    --max) max_parallel="${2:-}"; shift 2 ;;
    --target) target="${2:-}"; shift 2 ;;
    --interval) interval_seconds="${2:-}"; shift 2 ;;
    --lock-file) lock_file="${2:-}"; shift 2 ;;
    --once) once=1; shift ;;
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
if ! [[ "$interval_seconds" =~ ^[0-9]+$ ]] || [[ "$interval_seconds" -lt 1 ]]; then
  echo "--interval must be >= 1 seconds." >&2
  exit 1
fi
if [[ -z "$target" ]]; then
  target="$max_parallel"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
reconcile_script="${script_dir}/reconcile_ephemeral_pool.sh"
if [[ ! -x "$reconcile_script" ]]; then
  echo "missing executable: $reconcile_script" >&2
  exit 1
fi

run_reconcile() {
  local args=(
    "$reconcile_script"
    --repo "$repo"
    --min "$min_pool"
    --max "$max_parallel"
    --target "$target"
    --lock-file "$lock_file"
  )
  if [[ "$dry_run" -eq 1 ]]; then
    args+=(--dry-run)
  fi
  "${args[@]}"
}

if [[ "$once" -eq 1 ]]; then
  run_reconcile
  exit 0
fi

while true; do
  if ! run_reconcile; then
    echo "reconcile failed; retrying in ${interval_seconds}s" >&2
  fi

  if [[ "$dry_run" -eq 1 ]]; then
    break
  fi
  sleep "$interval_seconds"
done
