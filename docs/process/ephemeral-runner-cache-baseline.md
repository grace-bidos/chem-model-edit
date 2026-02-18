# Ephemeral Runner Cache Baseline

This document records reproducible CI baseline evidence for cache behavior and execution times.

Related issues:

- GitHub: #363
- Linear: GRA-100

## Scope

- Jobs: `web`, `api`, `contract`
- Workflow: `CI`
- Environment: self-hosted trusted runner lane

## Baseline sample (recorded)

- Date: 2026-02-18
- Run: <https://github.com/grace-bidos/chem-model-edit/actions/runs/22130319877>
- Evidence jobs:
  - web: <https://github.com/grace-bidos/chem-model-edit/actions/runs/22130319877/job/63968927164>
  - api: <https://github.com/grace-bidos/chem-model-edit/actions/runs/22130319877/job/63968927146>
  - contract: <https://github.com/grace-bidos/chem-model-edit/actions/runs/22130319877/job/63968927154>

Observed runner headers in logs:

- `Runner name: home-self-host*`
- `Machine name: grace-desktop`

## Baseline timings

All values are from GitHub Actions job metadata (`startedAt` -> `completedAt`) and key setup/test steps.

| Job | Total duration | Key cache/setup evidence |
| --- | --- | --- |
| `web` | 213s | `Restore pnpm store cache`: 102s, `Install dependencies`: 17s |
| `api` | 150s | `Setup uv`: 1s, `Install dependencies`: 1s, `Pytest`: 56s |
| `contract` | 157s | `Restore pnpm store cache`: 99s, `Install Python dependencies`: 1s, `Install Node dependencies`: 16s |

## Cache policy mapping

Current policy from `.github/workflows/ci.yml`:

- pnpm cache path:
  - `~/.pnpm-store`
- pnpm cache key:
  - `${{ runner.os }}-node22-pnpm-${{ hashFiles('pnpm-lock.yaml') }}`
- pnpm save policy:
  - save on `main` only when cache is missed
- uv cache:
  - enabled via `astral-sh/setup-uv` with restore/save for `apps/api/uv.lock`

## Reproduction steps

1. Select a completed CI run:

```bash
gh run list --workflow CI --limit 10
```

2. Inspect job durations and key steps:

```bash
gh run view <RUN_ID> --json jobs | jq -r '
  .jobs[]
  | select(.name=="web" or .name=="api" or .name=="contract")
  | .name as $n
  | "JOB \($n) total=\(((.completedAt|fromdateiso8601)-(.startedAt|fromdateiso8601)))s",
    (.steps[] | select(
      .name=="Restore pnpm store cache" or
      .name=="Install dependencies" or
      .name=="Install Python dependencies" or
      .name=="Install Node dependencies" or
      .name=="Setup uv" or
      .name=="Pytest"
    ) | "  STEP " + .name + " dur=\(((.completedAt|fromdateiso8601)-(.startedAt|fromdateiso8601)))s")
'
```

3. Confirm runner class (self-hosted vs hosted):

```bash
gh run view --job <JOB_ID> --log | head -n 12
```

## Update policy

Update this baseline when one of the following happens:

- cache key strategy changes in `.github/workflows/ci.yml`
- runner class changes (self-hosted <-> hosted) for trusted CI lane
- median job totals drift by more than 25% for two consecutive main-branch runs

When updating:

- append new run IDs and dates
- keep prior entries for trend visibility
