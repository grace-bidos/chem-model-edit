# Codex Agent Context (chem-model-edit)

## Project Summary

Web application for computational chemistry structure visualization and editing.
Primary capabilities include parallel structure views, partial structure transplant, compare/align, and export/share.

## Tech Stack

- Frontend: TanStack Start (SPA) + shadcn/ui + Mol\*
- Backend: FastAPI + ASE / pymatgen
- JS/TS: pnpm
- Python: uv (+ ruff, mypy, pytest)
- Tooling: Nx, Storybook, Chromatic
- Deploy: Cloudflare Workers (web), Cloud Run (api)

## Language and Communication

- Report to the user in Japanese.
- Public artifacts (GitHub PR/Issue/comments) must be written in English.
- Inter-agent communication (main agent <-> sub-agents, and sub-agent outputs) must be in English by default.
- Even when the initiating prompt is in Japanese, sub-agents should answer in English unless explicitly requested otherwise for a specific task.

## Mandatory Working Rules

- Always create and use a dedicated git worktree before starting implementation work.
- Keep `main` worktree for review/inspection only; do implementation in `.worktrees/<name>`.
- Do not develop directly on `main`.
- `uploads/` user-uploaded files may be moved/cleaned up if needed for the task.

## Planning and Delivery Model (Current)

- Source of truth for planning/progress: **Linear**.
- GitHub Issues are optional mirrors or external discussion threads.
- Delivery style: **Trunk-based + Stacked PR**.
- Branch format: `feature/identifier-title` (example: `feature/GRA-21-define-cutover-flags`).
- Planning hierarchy (fixed): `Project -> Parent Capability Issue -> Child Delivery Issue`.
- Parent capability issues track outcomes and acceptance criteria; do not use them as direct implementation containers.
- Child delivery issues are implementation units. One child issue maps to exactly one PR layer in a stack.

### Linear operation rules (fixed)

- Treat Linear Project/Issues as canonical. Do not use GitHub Issues as planning source-of-truth.
- Keep issue hierarchy to one level only (`Parent -> Child`). Do not create deeper nesting.
- Mandatory metadata for active parent and child issues (`Todo`/`In Progress`/`In Review`): `state`, `cycle`, `type:*`, and `size:*`.
- For delivery work in a project, child issues must have a parent capability issue. Top-level epics/capabilities are the only parentless issues in that project.
- Use `Backlog` only for parked work outside active execution. When an issue is pulled into the current cycle, move it to `Todo` or `In Progress`.
- Create follow-up child issues before closing an item when a PR is a partial/first slice.
- Do not close capability work by merge-only signal. Verify issue-level acceptance criteria first.

### Cycle assignment policy

- Put child delivery issues in cycles.
- Put parent capability issues in cycles by default to reflect current delivery velocity.
- Exception: keep a parent capability issue out of cycle only when it is intentionally long-running and not active this week.
- During cycle calibration phases, it is acceptable to intentionally over-pack current cycle scope to measure realistic weekly throughput.

### State transition policy

- Child issue lifecycle:
  - `Todo/In Progress/In Review`: active delivery states.
  - `Done`: allowed when PR is merged and child issue acceptance criteria are satisfied.
- Parent issue lifecycle:
  - `In Progress`: default when one or more child issues are active in the current cycle.
  - Move to `Done` only when required child issues are done and parent acceptance criteria are satisfied.
- If PR description includes wording like "first slice", "follow-up required", or deferred TODO boundaries, keep the parent open and create missing child issues immediately.

### Operational issue policy

- Do not create perpetual projects for ongoing environment/tooling chores.
- Default handling for operational chores: standalone team issue (no project), with `type:*` and `size:*`; assign to cycle only when actively worked.
- Promote operational work into a project only when it becomes a bounded, multi-issue initiative with explicit completion criteria.

### Work-item taxonomy

- Type axis: `Ask` / `Show` / `Ship`
  - Ask: spec/plan/architecture/DB/security-impacting decisions
  - Show: runtime behavior or feature changes
  - Ship: docs/tests/polish
- Size axis: `XS` / `S` / `M` / `L`
  - `L` must be split before merge queue entry

### Merge Queue and CI gate policy

- Queue required for runtime-impacting Ask and all Show changes.
- Queue optional for docs-only Ask/Ship items.
- Pre-merge required checks should stay fast: lint, typecheck, unit/smoke, policy checks.
- Merge-readiness gate (mandatory): required checks are green, unresolved review threads are zero, and head branch status is not `BEHIND` base.
- Contract-sensitive rule: when API contract changes, commit OpenAPI spec updates and regenerated client artifacts in the same PR.
- Heavy suites (if added later) should run post-merge.

## Tool Roles

- Linear: project/cycle/issue planning and status.
- Graphite (`gt`): stacked PR workflow (`create`, `submit`, `restack`, `sync`).
- GitHub CLI (`gh`): PR review/merge/check inspection.
- Git (`git`): low-level repository operations.
- Jujutsu (`jj`): optional local-history helper only; do not use as primary PR/push path.
- `scripts/gh/pr-autoloop.py`: lightweight PR loop helper for watch/blocker detection and optional auto-merge.

### PR Loop Operation (Recommended)

- Use `scripts/gh/pr-autoloop.py <PR_NUMBER> --watch --merge-when-ready --merge-method merge` to run the review/check/merge loop.
- Use `--resolve-outdated-threads` only for outdated threads; do not auto-resolve active review discussions.
- Prefer this script over ad-hoc manual polling when handling repeated CI/review feedback cycles.
- CodeRabbit trigger policy: request review once when PR is review-ready, then re-request only after substantive new commits and at least 15 minutes since the last request.
- If `--delete-branch` fails because the branch is attached to an active worktree, clean it up manually after removing the worktree.
- After a PR is merged, always sync local `main` immediately as part of the same operation set.
- Worktree-safe sync from any worktree:
  - `git fetch origin --prune`
  - `git -C "$(git rev-parse --show-toplevel)" merge --ff-only origin/main`

### Main vs Sub-agent execution model

- Assign sub-agents only child delivery issues, never parent capability issues.
- Main agent owns Linear planning/status/dependency management, parent acceptance checks, and final merge execution.
- Sub-agents own implementation, research, review-loop handling, and CI fixes for their assigned child issue and owned files.
- Inter-agent communication (including sub-agent outputs) is English by default.
- Parallel lane rule: one child issue maps to one lane (`issue -> branch -> worktree -> PR`); sub-agents must not edit outside their lane ownership.
- Lane conflict classification and default concurrency:
  - Low conflict lane (docs/specs/localized tests): up to 3 concurrent child lanes.
  - Medium conflict lane (same app area, clear file ownership): up to 2 concurrent child lanes.
  - High conflict lane (shared contracts/schemas/core runtime): 1 active lane; serialize merges.
- Sub-agent slot budget: up to 3 planned delivery lanes plus 1 reserved hotfix lane for `main` health recovery.
- Merge-readiness contract for sub-agent handoff: required checks are green, unresolved review threads are zero, head is not `BEHIND` base, and remaining risks/conflicts are listed.
- Conflict handoff rule: if cross-lane conflicts are detected, stop local conflict resolution and hand off to main agent with impacted files, conflict summary, and proposed resolution options.
- Sub-agents must not merge PRs directly. Main agent performs final merge and post-merge Linear updates.

## Repository Layout

- `apps/web`: TanStack Start frontend
- `apps/api`: FastAPI backend (`app/routers`, `services`, `tests`)
- `packages/api-client`: shared API client
- `packages/shared`: shared types
- `docs`: operational and process docs
- `specs`: spec/plan/task docs
- `scripts`: development helper scripts
- `samples/*`: QE/PDB/CIF samples
- `investigations`: research notes
- `ref-legacy`: legacy reference

## Development and Verification Commands

- Web:
  - `pnpm -C apps/web dev --port 3001`
  - `pnpm -C apps/web typecheck`
  - `pnpm -C apps/web storybook`
  - `pnpm -C apps/web chromatic`
- API:
  - `uv run uvicorn apps.api.main:app --reload --port 8000`
  - `uv run pytest`
  - `uv run mypy .`
- Cross-project:
  - `just style`
  - `just test`
  - `just typecheck`
  - `pnpm exec nx run-many -t lint,typecheck,test,knip`

## Domain/Implementation Notes

- Mol\* viewer center: `apps/web/src/components/molstar/MolstarViewer.tsx`
- Mol\* inputs: `pdbText` / `cifUrl` (single and multiple structures)
- Main API I/O assumption: Quantum ESPRESSO `.in`
- pnpm store policy: shared user-level store `~/.pnpm-store` (reuse across worktrees)

## Naming Rule for API Layers

- Distinguish layers by directory (`app/routers`, `services`).
- Keep the same domain filename across layers (for example, `routers/supercells.py` and `services/supercells.py`).
- Do not use singular/plural variation as a layer discriminator.

## Detailed Process References

Use these docs for full operational details instead of duplicating long instructions here:

- `docs/process/linear-stacked-pr-operating-model.md`
- `docs/process/coderabbit-parallel-playbook.md`
- `docs/process/merge-and-cleanup.md`
- `docs/process/gh-template-workflow.md`
