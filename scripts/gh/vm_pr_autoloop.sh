#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STACK_LOOP="${SCRIPT_DIR}/stack_lane_loop.py"

MAX_WAIT=5400
INTERVAL=20
MERGE_METHOD="merge"
GT_SYNC=1
RESOLVE_OUTDATED=0
CONFIRM_OUTDATED=0
DRY_RUN=0
PR_REF=""

usage() {
  cat <<'EOF'
Usage:
  scripts/gh/vm_pr_autoloop.sh <pr-number-or-url> [options]

Options:
  --max-wait <seconds>              Single watch cycle timeout (default: 5400)
  --interval <seconds>              Poll interval (default: 20)
  --merge-method <merge|squash|rebase>
  --no-gt-sync                      Skip gt sync
  --resolve-outdated-threads        Resolve only outdated unresolved threads
  --confirm-outdated-resolution     Required with --resolve-outdated-threads
  --dry-run                         Print command only
  -h, --help                        Show help
EOF
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "error: missing command: ${cmd}" >&2
    exit 1
  fi
}

while (($# > 0)); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --max-wait)
      MAX_WAIT="${2:-}"
      shift 2
      ;;
    --interval)
      INTERVAL="${2:-}"
      shift 2
      ;;
    --merge-method)
      MERGE_METHOD="${2:-}"
      shift 2
      ;;
    --no-gt-sync)
      GT_SYNC=0
      shift
      ;;
    --resolve-outdated-threads)
      RESOLVE_OUTDATED=1
      shift
      ;;
    --confirm-outdated-resolution)
      CONFIRM_OUTDATED=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -*)
      echo "error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      if [[ -n "${PR_REF}" ]]; then
        echo "error: PR reference already set: ${PR_REF}" >&2
        exit 1
      fi
      PR_REF="$1"
      shift
      ;;
  esac
done

if [[ -z "${PR_REF}" ]]; then
  echo "error: PR reference is required." >&2
  usage >&2
  exit 1
fi

if ! [[ "${MAX_WAIT}" =~ ^[0-9]+$ ]] || ((MAX_WAIT <= 0)); then
  echo "error: --max-wait must be a positive integer for single-cycle mode." >&2
  exit 1
fi

if ! [[ "${INTERVAL}" =~ ^[0-9]+$ ]] || ((INTERVAL < 5)); then
  echo "error: --interval must be an integer >= 5." >&2
  exit 1
fi

if [[ "${MERGE_METHOD}" != "merge" && "${MERGE_METHOD}" != "squash" && "${MERGE_METHOD}" != "rebase" ]]; then
  echo "error: --merge-method must be one of merge, squash, rebase." >&2
  exit 1
fi

if ((RESOLVE_OUTDATED == 1 && CONFIRM_OUTDATED == 0)); then
  echo "error: --resolve-outdated-threads requires --confirm-outdated-resolution." >&2
  exit 1
fi

if ((DRY_RUN == 0)); then
  require_cmd gh
  gh auth status >/dev/null
fi

if [[ ! -x "${STACK_LOOP}" ]]; then
  echo "error: executable not found: ${STACK_LOOP}" >&2
  exit 1
fi

cmd=(
  "${STACK_LOOP}" "${PR_REF}"
  --watch
  --merge-when-ready
  --merge-method "${MERGE_METHOD}"
  --interval "${INTERVAL}"
  --max-wait "${MAX_WAIT}"
)

if ((GT_SYNC == 1)); then
  cmd+=(--gt-sync)
fi

if ((RESOLVE_OUTDATED == 1)); then
  cmd+=(--resolve-outdated-threads)
fi

if ((DRY_RUN == 1)); then
  cmd+=(--dry-run)
fi

printf '$'
for part in "${cmd[@]}"; do
  printf ' %q' "${part}"
done
printf '\n'

exec "${cmd[@]}"
