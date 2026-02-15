# CodeRabbit Parallel Playbook

This playbook converts review wait time into delivery time by using stacked PRs and separate worktrees.

## Preconditions

- Every task has a child issue with clear acceptance criteria.
- Each child issue maps to one small PR.
- Worktree path is `.worktrees/<branch-name>`.
- Lane conflict classification and default concurrency:
  - Low conflict lane (docs/specs/localized tests): up to 3 concurrent child lanes.
  - Medium conflict lane (same app area, clear file ownership): up to 2 concurrent child lanes.
  - High conflict lane (shared contracts/schemas/core runtime): 1 active lane; serialize merges.
- Slot budget: up to 3 planned delivery lanes plus 1 reserved hotfix lane for `main` health recovery.

## Standard Flow

1. Open PR-A as draft from branch A.
2. Start PR-B immediately in another worktree from branch A (stacked).
3. When PR-A review arrives, fix PR-A first.
4. Merge PR-A with merge commit.
5. Restack/sync PR-B onto merged base (prefer `gt sync`; use rebase only if not using Graphite).
6. Repeat for PR-C, PR-D.

## Merge Gate (Mandatory)

Do not merge product-impacting PRs when required checks are red or incomplete.

- Product-impacting means API/auth/worker/runtime behavior changes.
- Required checks are all branch-protection required checks (including PR policy checks).
- Merge-ready means all of: required checks green, unresolved review threads = 0, and PR head is not `BEHIND` base.
- If CodeRabbit is marked required in the PR template, wait for review output (or explicitly document why unavailable) before merge.
- If API contract changes, commit OpenAPI updates and regenerated client artifacts in the same PR before requesting merge.
- If checks are flaky, re-run or fix within the same PR cycle. Do not defer known red checks to `main`.

## Start-Next-Task Gate (while waiting)

Start next task if all are true:

- Current PR is opened and CI is green or non-blocking.
- Remaining known fixes are low risk and local to current PR.
- Next task can progress at least 60% without the previous PR merged.

## CodeRabbit Trigger Policy (Rate-Limit Safe)

- Trigger once when a PR becomes review-ready.
- Re-trigger only after substantive new commits; do not re-trigger for metadata/comment-only updates.
- Keep at least 15 minutes between trigger comments.
- Batch fixes and push once per review cycle.
- If required review does not appear after 30 minutes, post one retry `@coderabbitai review`; if still blocked, document the blocker in PR.

## Post-Merge Recovery Protocol

If a product PR is merged with failing checks by mistake:

1. Open a recovery issue immediately.
2. Submit a hotfix PR to make `main` green first.
3. Backfill CodeRabbit review requests on the merged PR(s).
4. Create follow-up PRs for must-fix findings.
5. Update branch protection/rules if the merge was not blocked automatically.

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

- Link child issue and parent capability issue (or parent epic when applicable).
- Declare whether CodeRabbit is required or optional.
- If temporary behavior exists, document final target and follow-up PR.
