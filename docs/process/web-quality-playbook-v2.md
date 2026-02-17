# Web Quality Playbook v2

This playbook defines strict-but-pragmatic local quality operation for `apps/web`.

## Goals

- Keep CI required gates minimal and fast.
- Run deeper checks locally with clear command tiers.
- Expand scope and strictness without introducing unstable/noisy gates.

## Command tiers

- Quick: `just quality-quick`
  - lint, typecheck, unit tests.
- Standard: `just quality-standard`
  - quick + a11y + knip + fast-check + coverage.
- Deep: `just quality-deep`
  - standard + dependency graph + storybook build (+ chromatic when token exists) + mutation + playwright smoke.

## Current strictness baseline

- Vitest coverage thresholds (`apps/web/vite.config.ts`):
  - lines >= 90
  - functions >= 90
  - statements >= 90
  - branches >= 85
- Stryker thresholds (`apps/web/stryker.config.json`):
  - high: 88
  - low: 75
  - break: 68

## Promotion rubric (local -> CI candidate)

Promote a deep/local check to CI required gate only when all are true for at least 5 consecutive `main` runs:

1. Failure rate <= 2% and failures are actionable (non-flaky).
2. Median runtime is acceptable for PR feedback loop.
3. No persistent environment-coupled failures.

Until then, keep the check local-first.

## Tool-specific policy

- Playwright: local-only for now.
  - Use `pnpm -C apps/web test:e2e:smoke` by default.
  - If local browser/runtime is broken, treat as temporary local blocker and continue with non-E2E gates.
- Chromatic: optional in local deep flow.
  - Runs only when `CHROMATIC_PROJECT_TOKEN` is set.
- Dependency Cruiser:
  - Errors enforce architectural boundaries and cycles.
  - Orphan detection remains warning-level for visibility without hard-blocking.

## Mutation scope policy

- Keep mutation focused on high-value modules first.
- Exclude tests/stories/mocks from mutation targets.
- Expand scope only after improving survivor-heavy modules with deterministic tests.

## Reporting checklist (PR)

- Include executed command tier(s).
- Include coverage summary and mutation summary (score + top survivor clusters).
- Include skipped items and reasons (if any).
