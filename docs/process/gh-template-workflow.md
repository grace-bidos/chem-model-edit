# GitHub Template Workflow

Use body files instead of inline shell strings for issue/PR creation to avoid shell expansion bugs.

## Bootstrap PR Body Metadata

```bash
scripts/gh/bootstrap_pr_body.sh GRA-60 /tmp/pr_body.md --type Ship --size XS --queue Optional --stack Standalone --coderabbit Optional
```

This writes a template-friendly PR body with all metadata keys validated by `pr-policy`.

## Create Issue

```bash
scripts/gh/create_issue_from_template.sh "Title" /tmp/issue_body.md task
```

Use this only when creating an optional GitHub mirror/discussion issue.
Linear remains planning source-of-truth.

## Create PR

Preferred shortcut (bootstrap + create in one command):

```bash
just pr-open GRA-60 "ship: update runtime docs" type=Ship size=XS queue=Optional stack=Standalone coderabbit=Optional
```

For local-first flow, run this before opening the PR:

```bash
just pre-push-strict
```

Manual (body-file) flow:

```bash
scripts/gh/create_pr_from_template.sh <base-branch-or-main> <head-branch> "PR title" /tmp/pr_body.md --draft
```

## Check PR Readiness

```bash
scripts/gh/pr_readiness.py <pr-number-or-url>
```

The summary includes check status counts, unresolved review thread count, and
`mergeStateStatus` (including a `BEHIND` action hint).

## Optional: Tail Failed CI Logs

```bash
scripts/gh/ci_log_tail.sh <pr-number-or-url> --lines 200
```

Or target an explicit workflow run:

```bash
scripts/gh/ci_log_tail.sh --run-id <run-id> --lines 120
```

## Regenerate API Contract Artifacts

```bash
scripts/api/regenerate_contract.sh
```

This exports `openapi.json`, regenerates `packages/api-client/src/generated/schema.ts`,
and prints a focused `git diff` hint for review.

## Notes

- Keep public issue/PR text in English.
- Prefer markdown files committed in `specs/` or `/tmp/*.md` when drafting bodies.
- Ensure PR body includes required metadata fields checked by PR policy (Linear issue, type, size, queue policy, stack, and CodeRabbit policy).
- For stacked delivery, lane owners should use `scripts/gh/stack_lane_loop.py` and provide handoff status with `docs/process/subagent-lane-handoff-template.md`.
