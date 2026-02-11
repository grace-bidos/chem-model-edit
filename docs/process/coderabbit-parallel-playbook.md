# CodeRabbit Parallel Playbook

This playbook converts review wait time into delivery time by using stacked PRs and separate worktrees.

## Preconditions

- Every task has a child issue with clear acceptance criteria.
- Each child issue maps to one small PR.
- Worktree path is `.worktrees/<branch-name>`.

## Standard Flow

1. Open PR-A as draft from branch A.
2. Start PR-B immediately in another worktree from branch A (stacked).
3. When PR-A review arrives, fix PR-A first.
4. Merge PR-A with merge commit.
5. Rebase/update PR-B onto latest `main` (or onto merged A equivalent).
6. Repeat for PR-C, PR-D.

## Start-Next-Task Gate (while waiting)

Start next task if all are true:

- Current PR is opened and CI is green or non-blocking.
- Remaining known fixes are low risk and local to current PR.
- Next task can progress at least 60% without the previous PR merged.

## Rate Limit / Bot Delay Handling

- Wait at least 15 minutes before re-requesting review.
- Batch fixes and push once per review cycle.
- Use one follow-up comment summarizing all addressed items.
- If no review appears after 30 minutes, post `@coderabbitai review`.

## CodeRabbit Required vs Optional

Required:

- API contract changes
- AuthN/AuthZ logic changes
- Queue ownership or access control changes
- Worker execution behavior changes

Optional:

- Process docs
- Local helper scripts
- Template-only changes
- Non-runtime Justfile workflow changes

## PR Template Checklist

- Link child issue and parent epic.
- Declare whether CodeRabbit is required or optional.
- If temporary behavior exists, document final target and follow-up PR.
