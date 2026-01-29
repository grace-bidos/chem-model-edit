# Setup (local development)

This guide covers the local dev environment for `chem-model-edit`.
Quantum ESPRESSO (QE) runtime is **not** required here.
For the ZPE compute-plane (worker) setup, see `docs/zpe-worker-setup.md`.

## Prerequisites
- Git
- Node.js with Corepack enabled
- pnpm 10.27.0
- Python >= 3.13 (use `uv`)
- Optional: `just` (recommended for `just dev`)

If you want `just`:
- With Rust/Cargo installed: `cargo install just`
- Or via your OS package manager

## Quick start (script)
```bash
./scripts/setup-dev.sh
```
Notes:
- The script runs `pnpm approve-builds` and asks you to approve `@parcel/watcher`.
- `pnpm approve-builds` will update `pnpm-workspace.yaml`.
- Set `SETUP_SKIP_APPROVE=1` to skip the approval step.

## Manual setup
```bash
# Node + pnpm
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

# Approve build scripts (select @parcel/watcher)
pnpm approve-builds

# Python + API deps
uv python install 3.13
cd apps/api
uv sync
```

## Run
```bash
# Web
pnpm -C apps/web dev --port 3000

# API
cd apps/api
uv run uvicorn main:app --reload --port 8000
```

Run both together (requires `just`):
```bash
just deps
just dev
```

## Quality checks
```bash
pnpm -C apps/web typecheck
pnpm -C apps/web lint
pnpm -C apps/web test

cd apps/api
uv run pytest
uv run mypy .
```

Just recipes:
```bash
just style      # web: prettier + eslint, api: ruff
just test       # web: vitest, api: pytest
just typecheck  # web: tsc, api: mypy
```

## Troubleshooting
- In sandboxed environments (for example, Codex CLI), run `uv sync` with elevation.
  Without elevation, `uv` can fail with `EXDEV` (cross-device rename).
