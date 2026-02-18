# Merge And Cleanup Workflow

Use this flow to avoid branch deletion failures when a branch is still attached to a worktree.
Prefer the PR auto-loop for routine operations; use manual steps as fallback.

## Rules

- Standard split: use `gt` only for stack operations (`sync/restack/submit`) and `gh` scripts for readiness and merge loop.
- Prefer `scripts/gh/stack_lane_loop.py` for lane-owner routine flow (`--gt-sync` + readiness loop).
- `scripts/gh/pr-autoloop.py` remains the direct watcher and merge executor.
- Lane owner (sub-agent) runs the watch loop and CI/review handling for that lane.
- Main agent should not continuously poll PR/CI for all lanes; it executes milestone-based checks and final merge.
- In a new worktree, run preflight checks before watch mode: `for tool in node corepack pnpm uv; do command -v "$tool" >/dev/null || { echo "missing: $tool"; exit 1; }; done`, `node -v`, `pnpm -v`, `uv --version`.
- Merge-ready is based on required checks + review threads + branch-behind status. Non-required pending checks do not block merge.
- If `--delete-branch` fails because branch is attached to a worktree, use a two-step close: merge first, then cleanup worktree/branch.
- Always sync local `main` immediately after merge.
- On merge/reassignment/cancel, close lane explicitly: post final status, clean worktree/branch, and free slot inventory.

## Commands

```bash
# 0) Preferred lane-owner entrypoint (gt sync + readiness/watch/merge)
scripts/gh/stack_lane_loop.py <PR_NUMBER> --gt-sync --watch --merge-when-ready --merge-method merge

# 1) Optional one-shot readiness summary
scripts/gh/pr_readiness.py <PR_NUMBER_OR_URL>

# 2) Direct watcher/merge command (same core loop used by helper entrypoint)
scripts/gh/pr-autoloop.py <PR_NUMBER> --watch --merge-when-ready --merge-method merge

# 3) Optional manual fallback from main worktree
scripts/git/merge_pr.sh <pr-number-or-url>

# 4) Cleanup worktree and local branch (if merge used no --delete-branch)
scripts/git/cleanup_worktree.sh .worktrees/<name> <branch>

# 5) Sync local main
git fetch origin --prune
git -C "$(git rev-parse --show-toplevel)" merge --ff-only origin/main
```

If checks are failing and you need quick context, tail failed-step logs:

```bash
scripts/gh/ci_log_tail.sh <PR_NUMBER_OR_URL> --lines 200
```

## Dry-run

```bash
scripts/git/cleanup_worktree.sh --dry-run .worktrees/<name> <branch>
```
