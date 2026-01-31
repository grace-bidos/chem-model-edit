set shell := ["bash", "-cu"]
# Note: In sandboxed environments (e.g. Codex CLI), run these recipes with elevation.
web_port := "3001"
api_port := "8000"

@default:
  just --list

setup:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/setup-dev.sh

api-sync:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv sync
  popd >/dev/null

deps:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm install
  just api-sync

web:
  #!/usr/bin/env bash
  set -euo pipefail
  port="${WEB_PORT:-{{web_port}}}"
  initial_port="$port"
  port="$(./scripts/find_free_port.sh "$port")"
  if [[ "$port" != "$initial_port" ]]; then
    echo "WEB_PORT $initial_port は使用中のため $port を利用します"
  fi
  pnpm -C apps/web dev --port "$port"

api:
  #!/usr/bin/env bash
  set -euo pipefail
  port="${API_PORT:-{{api_port}}}"
  initial_port="$port"
  port="$(./scripts/find_free_port.sh "$port")"
  if [[ "$port" != "$initial_port" ]]; then
    echo "API_PORT $initial_port は使用中のため $port を利用します"
  fi
  cd apps/api
  uv run uvicorn main:app --reload --port "$port"

dev:
  #!/usr/bin/env bash
  set -euo pipefail
  trap 'kill 0' INT TERM EXIT
  api_port="${API_PORT:-{{api_port}}}"
  web_port="${WEB_PORT:-{{web_port}}}"
  initial_api_port="$api_port"
  initial_web_port="$web_port"
  api_port="$(./scripts/find_free_port.sh "$api_port")"
  web_port="$(./scripts/find_free_port.sh "$web_port")"
  if [[ "$api_port" != "$initial_api_port" ]]; then
    echo "API_PORT $initial_api_port は使用中のため $api_port を利用します"
  fi
  if [[ "$web_port" != "$initial_web_port" ]]; then
    echo "WEB_PORT $initial_web_port は使用中のため $web_port を利用します"
  fi
  pushd apps/api >/dev/null
  uv run uvicorn main:app --reload --port "$api_port" &
  popd >/dev/null
  pnpm -C apps/web dev --port "$web_port" &
  wait

zpe-worker:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/run-zpe-worker.sh

zpe-http-worker:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/run-zpe-http-worker.sh

api-test:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run pytest
  popd >/dev/null

api-ruff:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run ruff check .
  popd >/dev/null

api-mypy:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run mypy .
  popd >/dev/null

web-format:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web format

web-lint:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web lint

web-typecheck:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web typecheck

web-test:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test

web-test-ui:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:ui

storybook:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web storybook

storybook-build:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web build-storybook

chromatic:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web chromatic

nx:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx

nx-graph:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx graph

nx-storybook:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run web:storybook

nx-chromatic:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run web:chromatic

nx-lint:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t lint

nx-typecheck:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t typecheck

nx-test:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t test

nx-build:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t build

nx-knip:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t knip

style:
  #!/usr/bin/env bash
  set -euo pipefail
  just web-format
  just web-lint
  just api-ruff

test:
  #!/usr/bin/env bash
  set -euo pipefail
  just web-test
  just api-test

typecheck:
  #!/usr/bin/env bash
  set -euo pipefail
  just web-typecheck
  just api-mypy

ci:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t lint,typecheck,test,knip
