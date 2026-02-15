# Codex Agent Context (chem-model-edit)

## Project Summary
Web application for computational chemistry structure visualization and editing.
Primary capabilities include parallel structure views, partial structure transplant, compare/align, and export/share.

## Tech Stack
- Frontend: TanStack Start (SPA) + shadcn/ui + Mol*
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
- One implementation issue maps to one PR layer in a stack.

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
- If `--delete-branch` fails because the branch is attached to an active worktree, clean it up manually after removing the worktree.
- After a PR is merged, always sync local `main` immediately as part of the same operation set.
- Worktree-safe sync from any worktree:
  - `git fetch origin --prune`
  - `git -C "$(git rev-parse --show-toplevel)" merge --ff-only origin/main`

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
- Mol* viewer center: `apps/web/src/components/molstar/MolstarViewer.tsx`
- Mol* inputs: `pdbText` / `cifUrl` (single and multiple structures)
- Main API I/O assumption: Quantum ESPRESSO `.in`
- pnpm store policy: `.pnpm-store`

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
