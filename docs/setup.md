# Setup (local development)

This guide covers the local dev environment for `chem-model-edit`.
Quantum ESPRESSO (QE) runtime is **not** required here.
For the ZPE compute-plane (worker) setup, see `docs/zpe-worker-setup.md`.

## Prerequisites
- Git
- Node.js with Corepack enabled
- pnpm 10.27.0
- Python >= 3.13 (use `uv`)
- Optional: `just` (`just dev` for Web+API together)

`just` is recommended (but optional) for daily development commands like `just dev`, `just test`, and `just typecheck`.

If you want `just`:
- Via Cargo (requires Rust/Cargo environment): `cargo install just`
- On Linux via snap: `sudo snap install --classic just`
- Or via your OS package manager

## Recommended path
`./scripts/setup-dev.sh` is the recommended entrypoint.

```bash
./scripts/setup-dev.sh
```

What this setup does:
- Runs `corepack enable` + `corepack prepare pnpm@10.27.0 --activate`
- Runs `pnpm install`
- Runs `pnpm approve-builds` (approve `@parcel/watcher`)
- Runs `uv python install 3.13` and `uv sync` in `apps/api`

Notes:
- `pnpm approve-builds` updates `pnpm-workspace.yaml`.
- Set `SETUP_SKIP_APPROVE=1` to skip `pnpm approve-builds`.
- Set `SETUP_INSTALL_JUST=1` to auto-install `just` via `cargo` (requires Rust/Cargo, and only runs if `cargo` exists).

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
### Start both Web and API (recommended)
```bash
just dev
```

`just dev` requires `just` to be installed.

Port behavior for `just dev`:
- Default values are `WEB_PORT=3001` and `API_PORT=8000`.
- If a port is already in use, `just dev` automatically selects a free port and prints it.

Override ports explicitly:
```bash
WEB_PORT=4000 API_PORT=9000 just dev
```

### Start components separately (alternative)
```bash
# Web
pnpm -C apps/web dev

# API
cd apps/api
uv run uvicorn main:app --reload --port 8000
```
When Web is started directly with Vite, the default port is `3000`.

## Quality checks
```bash
just style
just typecheck
just test
just ci
```

Recipe details: `just --list` and `Justfile`.

## Troubleshooting
- In sandboxed environments (for example, Codex CLI), run `uv sync` with elevation.
  Without elevation, `uv` can fail with `EXDEV` (cross-device rename).
