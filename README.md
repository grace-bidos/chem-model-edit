# chem-model-edit

English | [Japanese](./README.ja.md)

[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/grace-bidos/chem-model-edit)

Web app for editing and visualizing chemical structures (Quantum ESPRESSO `.in`).

Goal: build a browser workflow to view, edit, compare, and share QE `.in` structures.
Detailed scope and acceptance criteria (archived reference): `specs/archive/001-chem-model-webapp/spec.md`.
Current planning source-of-truth: Linear. Specs directory guide: `specs/README.md`.

## Quickstart

Run from repository root.

```bash
git clone git@github.com:grace-bidos/chem-model-edit.git
cd chem-model-edit

./scripts/setup-dev.sh
```

`just` is recommended for daily development (`just dev`, `just test`, `just typecheck`), but optional.

Install `just` manually:

```bash
# Requires Rust/Cargo environment
cargo install just

# Linux (snap)
sudo snap install --classic just
```

No `just`? start components separately:

```bash
pnpm -C apps/web dev
cd apps/api && uv run uvicorn main:app --reload --port 8000
```

If you want to auto-install `just` via `cargo` during setup:

```bash
SETUP_INSTALL_JUST=1 ./scripts/setup-dev.sh
just dev
```

- `./scripts/setup-dev.sh` installs pnpm deps, approves build scripts, and runs `uv sync` for API.
- `just` is optional.
- `just dev` starts Web + API together.
  - Default ports via `just`: `WEB_PORT=3001`, `API_PORT=8000`.
  - If a port is in use, `just dev` auto-selects a free port and prints it.

Full setup details: `docs/setup.md`  
ZPE worker setup: `docs/zpe-worker-setup.md`

## Run Components Separately (alternative)

All commands below run from repository root unless noted.

### Web only

```bash
pnpm -C apps/web dev
```

Default port is `3000` for direct Vite dev.

### API only

```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload --port 8000
```

### Override ports for `just dev`

```bash
WEB_PORT=4000 API_PORT=9000 just dev
```

## Quality checks

```bash
just style      # web: prettier + eslint, api: ruff
just typecheck  # web: tsc, api: mypy
just test       # web: vitest, api: pytest
just ci         # nx run-many -t lint,typecheck,test,knip
just quality-quick     # local quick gate: lint + typecheck + unit
just ci-pr-quick       # PR相当の quick gate: web/api/contract drift
just git-hooks-install # pre-push strict hook を有効化
just quality-standard  # local standard gate: + a11y + knip + fastcheck + coverage
just quality-deep      # local deep gate: + depcruise + storybook/chromatic + stryker + playwright smoke
```

Detailed policy: `docs/process/web-quality-playbook-v2.md`
CI/CD運用: `docs/process/local-first-ci-cd.md`

## Repository layout

- `apps/web`: TanStack Start web app (frontend)
- `apps/api`: FastAPI backend
- `packages/api-client`: shared API client package
- `packages/shared`: shared types
- `docs`: setup and operations docs
- `specs`: spec/plan/task documents
- `samples/qe-in`: sample QE `.in` inputs

## Notes (sandboxed environments)

- In sandboxed environments (for example, Codex CLI), `uv sync` or `just` may require elevated execution due to filesystem constraints (`EXDEV`).
