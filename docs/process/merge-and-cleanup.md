# Merge And Cleanup Workflow

Use this flow to avoid branch deletion failures when a branch is still attached to a worktree.
Prefer the PR auto-loop for routine operations; use manual steps as fallback.

## Rules

- Prefer `scripts/gh/pr-autoloop.py` with `--merge-when-ready`.
- If `--delete-branch` fails because branch is attached to a worktree, remove worktree/branch manually afterward.
- Always sync local `main` immediately after merge.

## Commands

```bash
# 1) Preferred: automated review/check/merge loop
scripts/gh/pr-autoloop.py <PR_NUMBER> --watch --merge-when-ready --merge-method merge

# 2) Optional manual fallback from main worktree
scripts/git/merge_pr.sh <pr-number-or-url>

# 3) Cleanup worktree and local branch
scripts/git/cleanup_worktree.sh .worktrees/<name> <branch>

# 4) Sync local main
git fetch origin --prune
git -C "$(git rev-parse --show-toplevel)" merge --ff-only origin/main
```

## Dry-run

```bash
scripts/git/cleanup_worktree.sh --dry-run .worktrees/<name> <branch>
```
