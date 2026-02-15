#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
usage:
  scripts/gh/ci_log_tail.sh <PR_NUMBER_OR_URL> [--lines N]
  scripts/gh/ci_log_tail.sh --run-id <RUN_ID> [--lines N]

Tail failed-step logs from a workflow run.
If PR is provided, the latest non-success run on that PR head branch is selected.
USAGE
}

lines=120
run_id=""
pr_ref=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --run-id)
      run_id="${2:-}"
      shift 2
      ;;
    --lines)
      lines="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      if [[ -n "$pr_ref" ]]; then
        echo "error: unexpected argument: $1" >&2
        usage >&2
        exit 1
      fi
      pr_ref="$1"
      shift
      ;;
  esac
done

if [[ -z "$run_id" && -z "$pr_ref" ]]; then
  usage >&2
  exit 1
fi

if [[ -n "$run_id" && ! "$run_id" =~ ^[0-9]+$ ]]; then
  echo "error: --run-id must be numeric" >&2
  exit 1
fi

if [[ ! "$lines" =~ ^[0-9]+$ ]]; then
  echo "error: --lines must be numeric" >&2
  exit 1
fi

if [[ -z "$run_id" ]]; then
  head_branch="$(gh pr view "$pr_ref" --json headRefName --jq .headRefName)"
  run_id="$(gh run list --branch "$head_branch" --limit 30 --json databaseId,status,conclusion --jq '
    map(select((.status != "completed") or (.conclusion != "success")))
    | .[0].databaseId // ""
  ')"

  if [[ -z "$run_id" || "$run_id" == "null" ]]; then
    echo "No non-success workflow run found for branch: $head_branch" >&2
    exit 1
  fi
fi

echo "Showing failed-step logs for run id: $run_id (last $lines lines)"
gh run view "$run_id" --log-failed | tail -n "$lines"
