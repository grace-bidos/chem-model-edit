#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <pr-number-or-url>" >&2
  exit 1
fi

pr_ref="$1"
current_branch="$(git rev-parse --abbrev-ref HEAD)"

if [[ "$current_branch" != "main" ]]; then
  echo "error: run merge_pr.sh from main worktree only (current: $current_branch)" >&2
  echo "hint: switch to your main worktree and run again" >&2
  exit 1
fi

if ! git diff --quiet || ! git diff --cached --quiet; then
  echo "error: working tree is not clean. commit/stash before merging." >&2
  exit 1
fi

echo "Merging PR with merge commit: $pr_ref"
gh pr merge "$pr_ref" --merge
echo "Merge complete."
echo "Next: run scripts/git/cleanup_worktree.sh <worktree-path> <branch>"
