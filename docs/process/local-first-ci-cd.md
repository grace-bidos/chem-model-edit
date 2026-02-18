# Local-First CI/CD Operations

## Objective

Keep feedback loops fast by running quick checks locally, while keeping required CI lanes lean and deterministic.

## PR required checks (quick lanes)

- `web`: lint, typecheck, unit tests
- `api`: ruff, mypy, pytest + diff gates
- `contract`: OpenAPI drift + generated client drift

Heavy checks are not required in normal PR flow.

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
