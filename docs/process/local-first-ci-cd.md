# Local-First CI/CD Operations

## Objective

Keep feedback loops fast by running quick checks locally, while keeping required CI lanes lean and deterministic.

## PR required checks (quick lanes)

- `web`: lint, typecheck, unit tests
- `api`: ruff, mypy, pytest + diff gates
- `contract`: OpenAPI drift + generated client drift

Heavy checks are not required in normal PR flow.

Merge-ready is defined by required checks only:

- required checks are green
- unresolved review threads are zero
- PR head is not `BEHIND` base

`CodeRabbit` pending status is not a merge blocker unless it is explicitly configured as a required check.

## PR creation (metadata-safe)

Use the template-driven helper to avoid `pr-metadata` failures:

```bash
just pr-open GRA-208 "docs: retire cloud run legacy references"
```

Optional flags:

```bash
just pr-open GRA-208 "docs: retire cloud run legacy references" type=Ship size=S queue=Optional stack=Standalone coderabbit=Optional draft=true
```

Avoid direct `gh pr create` unless you pass a validated body file.

## Full checks (on demand)

Use the `Full CI (On Demand)` workflow for heavy suites.

Triggers:

- `workflow_dispatch`
- `schedule`
- PR labeled `ci:full`

Includes:

- web a11y/fastcheck
- api pyright/bandit/pip-audit/schemathesis broad
- optional projection smoke

## Cloudflare Workers deploy policy

Production deploy is label-gated.

- Add label `release:workers` to PR
- Merge PR into `main`
- Workflow `Deploy Web to Cloudflare Workers` runs automatically

Manual deploy remains available via `workflow_dispatch`.

## Docs-only PR behavior

On pull requests, quick lanes are routed by changed paths.
Docs-only changes should skip `web/api/contract` lanes unless CI workflow files affecting quick lanes changed.
On `push` to `main` and `merge_group`, quick lanes still run by design.

## Local strict gate

Install local hooks once:

```bash
just git-hooks-install
```

Then every push runs:

```bash
just pre-push-strict
```

Quick bypass for emergency cases only:

```bash
SKIP_PREPUSH_STRICT=1 git push
```

## Stacked lane quick flow

For each child-issue lane:

```bash
just pre-push-strict
just pr-open GRA-XXX "ship: concise PR title"
scripts/gh/stack_lane_loop.py <PR_NUMBER> --gt-sync --watch --merge-when-ready --merge-method merge
```

Before handoff to the main agent, fill:

```bash
just lane-handoff-template
```

## Weekly operation KPIs (2-sprint calibration)

Capture weekly values in Linear comments or a cycle note:

- median time-to-green for required checks (`web`, `api`, `contract`)
- PR lead time from open to merge
- CI rerun count per PR
- percentage of PRs merged without manual polling outside lane scripts

Use these KPIs to tune lane concurrency and identify slow checks without expanding required CI scope.
