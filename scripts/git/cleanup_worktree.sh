#!/usr/bin/env bash
set -euo pipefail

dry_run=0
if [[ "${1:-}" == "--dry-run" ]]; then
  dry_run=1
  shift
fi

if [[ $# -lt 2 ]]; then
  echo "usage: $0 [--dry-run] <worktree-path> <branch>" >&2
  exit 1
fi

worktree_path="$1"
branch="$2"

run() {
  if [[ $dry_run -eq 1 ]]; then
    echo "[dry-run] $*"
  else
    "$@"
  fi
}

echo "Cleaning worktree and branch"
run git worktree remove "$worktree_path"
run git worktree prune
run git branch -d "$branch"
echo "Cleanup complete: $worktree_path / $branch"
