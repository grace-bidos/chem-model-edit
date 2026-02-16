# VM Stack PR Auto-loop Ops Runbook

This runbook standardizes safe lane-owner operation for VM stack PR loops using `scripts/gh/pr-autoloop.py`.

## Scope

- Intended for child delivery issue lanes in the VM stack.
- Single-cycle watch mode only (finite timeout, no endless polling loops).
- Goal: merge only when merge-readiness gates are satisfied.

## Merge-ready criteria

Treat a PR as merge-ready only when all of the following are true:

- required checks are green
- unresolved review threads are `0`
- head branch is not `BEHIND` base
- PR state is `OPEN` and not draft
- `mergeStateStatus` is `CLEAN`

`scripts/gh/stack_lane_loop.py` runs `pr_readiness.py` first, then passes control to `pr-autoloop.py` watch mode.

## Outdated-thread handling policy

Default policy: do not auto-resolve any review thread.

Use `--resolve-outdated-threads` only when all of the following are true:

- every targeted unresolved thread is clearly marked outdated
- no active reviewer conversation is still open
- you verified no hidden follow-up work is required

Never auto-resolve non-outdated unresolved threads.

## Recommended command

Use the VM-safe wrapper (single-cycle defaults included):

```bash
scripts/gh/vm_pr_autoloop.sh <PR_NUMBER_OR_URL>
```

Defaults:

- `--watch --merge-when-ready`
- `--gt-sync`
- `--merge-method merge`
- `--interval 20`
- `--max-wait 5400` (90 minutes, single cycle)

## Command examples

Standard lane flow:

```bash
scripts/gh/vm_pr_autoloop.sh 123
```

Longer single cycle (2 hours):

```bash
scripts/gh/vm_pr_autoloop.sh 123 --max-wait 7200
```

Skip `gt sync` when stack/base is already aligned:

```bash
scripts/gh/vm_pr_autoloop.sh 123 --no-gt-sync
```

Outdated-only thread cleanup (explicit opt-in):

```bash
scripts/gh/vm_pr_autoloop.sh 123 \
  --resolve-outdated-threads \
  --confirm-outdated-resolution
```

Dry-run preview:

```bash
scripts/gh/vm_pr_autoloop.sh 123 --dry-run
```

## Guardrails

- Do not run infinite watch loops for VM lanes (`--max-wait` must stay positive).
- Do not use auto-resolve for active reviewer conversations.
- If merge is blocked by checks/reviews, stop and fix blockers first; rerun loop after updates.
- After merge, sync local `main` immediately per `docs/process/merge-and-cleanup.md`.
