# Merge And Cleanup Workflow

Use this flow to avoid branch deletion failures when a branch is still attached to a worktree.

## Rules

- Run PR merge from the `main` worktree only.
- Do not use `--delete-branch` in `gh pr merge`.
- Cleanup local worktree and branch after merge in a separate step.

## Commands

```bash
# 1) Merge the PR from main worktree
scripts/git/merge_pr.sh <pr-number-or-url>

# 2) Cleanup worktree and local branch
scripts/git/cleanup_worktree.sh .worktrees/<name> <branch>
```

## Dry-run

```bash
scripts/git/cleanup_worktree.sh --dry-run .worktrees/<name> <branch>
```
