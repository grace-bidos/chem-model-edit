#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "==> Node + pnpm"
if ! command -v corepack >/dev/null 2>&1; then
  echo "corepack not found. Install Node.js (with Corepack) first." >&2
  exit 1
fi
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

if [[ "${SETUP_SKIP_APPROVE:-}" != "1" ]]; then
  if grep -q "@parcel/watcher" pnpm-workspace.yaml; then
    echo "pnpm approve-builds already applied for @parcel/watcher."
  else
    echo "==> Approve build scripts (select @parcel/watcher)"
    pnpm approve-builds
  fi
fi

if command -v uv >/dev/null 2>&1; then
  echo "==> Python + API deps (uv)"
  uv python install 3.13
  (cd apps/api && uv sync)
else
  echo "uv not found. Install uv to set up Python deps." >&2
fi

if ! command -v just >/dev/null 2>&1; then
  if [[ "${SETUP_INSTALL_JUST:-}" == "1" ]]; then
    echo "==> Optional: install just"
    if command -v cargo >/dev/null 2>&1; then
      cargo install just
    else
      echo "cargo not found. Skip auto install of just." >&2
      echo "Install just manually if needed (e.g. cargo install just)." >&2
    fi
  else
    echo "just not found (optional)."
    echo "To auto-install via cargo: SETUP_INSTALL_JUST=1 ./scripts/setup-dev.sh"
  fi
fi

echo "Setup done."
