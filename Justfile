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

api-openapi-export:
  #!/usr/bin/env bash
  set -euo pipefail
  PYTHONPATH=apps/api uv run --project apps/api python apps/api/scripts/export_openapi.py

api-openapi-generate:
  #!/usr/bin/env bash
  set -euo pipefail
  just api-openapi-export
  mkdir -p docs/api
  rm -rf docs/api/openapi-html docs/api/openapi-markdown
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD:/local" \
    openapitools/openapi-generator-cli:v7.16.0 generate \
    -i /local/packages/api-client/openapi/openapi.json \
    -g html2 \
    -o /local/docs/api/openapi-html
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    -v "$PWD:/local" \
    openapitools/openapi-generator-cli:v7.16.0 generate \
    -i /local/packages/api-client/openapi/openapi.json \
    -g markdown \
    -o /local/docs/api/openapi-markdown
  echo "generated: docs/api/openapi-html/index.html"
  echo "generated: docs/api/openapi-markdown/README.md"

api-pyreverse:
  #!/usr/bin/env bash
  set -euo pipefail
  mkdir -p docs/graphs
  pushd apps/api >/dev/null
  uv run pyreverse -o dot -p api app services
  popd >/dev/null
  mv -f apps/api/classes_api.dot docs/graphs/api-classes.dot
  mv -f apps/api/packages_api.dot docs/graphs/api-packages.dot
  if command -v dot >/dev/null 2>&1; then
    dot -Tsvg docs/graphs/api-classes.dot -o docs/graphs/api-classes.svg
    dot -Tsvg docs/graphs/api-packages.dot -o docs/graphs/api-packages.svg
    echo "generated: docs/graphs/api-classes.svg"
    echo "generated: docs/graphs/api-packages.svg"
  else
    echo "skip svg export: graphviz 'dot' not found"
  fi
  echo "generated: docs/graphs/api-classes.dot"
  echo "generated: docs/graphs/api-packages.dot"

api-cov:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run pytest --cov=app --cov=services --cov-report=term-missing --cov-report=html:htmlcov --cov-report=xml:coverage.xml
  popd >/dev/null
  echo "generated: apps/api/htmlcov/index.html"
  echo "generated: apps/api/coverage.xml"

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

story-gen name dir export_style='named' export_name='' title_prefix='':
  #!/usr/bin/env bash
  set -euo pipefail
  name="{{name}}"
  dir="{{dir}}"
  export_style="{{export_style}}"
  export_name="{{export_name}}"
  title_prefix="{{title_prefix}}"
  args=(story new --name "$name" --dir "$dir" --exportStyle "$export_style")
  if [[ -n "$export_name" ]]; then
    args+=(--exportName "$export_name")
  fi
  if [[ -n "$title_prefix" ]]; then
    args+=(--titlePrefix "$title_prefix")
  fi
  pnpm exec hygen "${args[@]}"

story-gen-ui-batch:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/storybook/generate-ui-stories.sh

story-gen-editor-v2-batch:
  #!/usr/bin/env bash
  set -euo pipefail
  ./scripts/storybook/generate-editor-v2-stories.sh

story-catalog-check:
  #!/usr/bin/env bash
  set -euo pipefail
  url="http://127.0.0.1:6006"
  cleanup() {
    if [[ -n "${storybook_pid:-}" ]]; then
      kill "$storybook_pid" 2>/dev/null || true
      wait "$storybook_pid" 2>/dev/null || true
    fi
    pids="$(lsof -ti tcp:6006 -sTCP:LISTEN 2>/dev/null || true)"
    if [[ -n "$pids" ]]; then
      kill $pids 2>/dev/null || true
    fi
  }
  trap cleanup EXIT
  pnpm -C apps/web storybook --port 6006 >/tmp/storybook.log 2>&1 &
  storybook_pid=$!
  for i in {1..60}; do
    if curl -sf "$url" >/dev/null; then
      echo "Storybook is ready at $url"
      exit 0
    fi
    sleep 1
  done
  echo "Storybook did not become ready within timeout." >&2
  echo "--- storybook log ---" >&2
  cat /tmp/storybook.log >&2
  exit 1

nx:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx

nx-graph:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx graph

graph-nx:
  #!/usr/bin/env bash
  set -euo pipefail
  mkdir -p docs/graphs
  NX_DAEMON=false NX_ISOLATE_PLUGINS=false pnpm exec nx graph --file docs/graphs/nx-project-graph.html --open=false
  echo "generated: docs/graphs/nx-project-graph.html"

graph-web-deps:
  #!/usr/bin/env bash
  set -euo pipefail
  mkdir -p docs/graphs
  pnpm exec depcruise --config .dependency-cruiser.cjs --include-only "^apps/web/src" --output-type dot --output-to docs/graphs/web-dependency-graph.dot apps/web/src
  pnpm exec depcruise --config .dependency-cruiser.cjs --include-only "^apps/web/src" --output-type mermaid --output-to docs/graphs/web-dependency-graph.mmd apps/web/src
  pnpm exec depcruise --config .dependency-cruiser.cjs --include-only "^apps/web/src" --output-type json --output-to docs/graphs/web-dependency-graph.json apps/web/src
  if command -v dot >/dev/null 2>&1; then
    dot -Tsvg docs/graphs/web-dependency-graph.dot -o docs/graphs/web-dependency-graph.svg
    echo "generated: docs/graphs/web-dependency-graph.svg"
  else
    echo "skip svg export: graphviz 'dot' not found"
  fi
  echo "generated: docs/graphs/web-dependency-graph.dot"
  echo "generated: docs/graphs/web-dependency-graph.mmd"
  echo "generated: docs/graphs/web-dependency-graph.json"

graph-editor-v2-madge:
  #!/usr/bin/env bash
  set -euo pipefail
  mkdir -p docs/graphs
  madge_bin="./node_modules/.bin/madge"
  src_dir="apps/web/src/features/editor-v2/components"
  exclude_pattern="(\\.stories\\.|\\.test\\.|\\.a11y\\.test\\.|\\.fastcheck\\.test\\.)"
  "$madge_bin" "$src_dir" \
    --extensions ts,tsx \
    --ts-config apps/web/tsconfig.json \
    --exclude "$exclude_pattern" \
    --json > docs/graphs/editor-v2-madge.json
  "$madge_bin" "$src_dir" \
    --extensions ts,tsx \
    --ts-config apps/web/tsconfig.json \
    --exclude "$exclude_pattern" \
    --dot > docs/graphs/editor-v2-madge.dot
  "$madge_bin" "$src_dir" \
    --extensions ts,tsx \
    --ts-config apps/web/tsconfig.json \
    --exclude "$exclude_pattern" \
    --circular \
    --json > docs/graphs/editor-v2-madge.circular.json
  if command -v dot >/dev/null 2>&1; then
    "$madge_bin" "$src_dir" \
      --extensions ts,tsx \
      --ts-config apps/web/tsconfig.json \
      --exclude "$exclude_pattern" \
      --image docs/graphs/editor-v2-madge.svg
    echo "generated: docs/graphs/editor-v2-madge.svg"
  else
    echo "skip svg export: graphviz 'dot' not found"
  fi
  echo "generated: docs/graphs/editor-v2-madge.json"
  echo "generated: docs/graphs/editor-v2-madge.dot"
  echo "generated: docs/graphs/editor-v2-madge.circular.json"

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

git-main-sync:
  #!/usr/bin/env bash
  set -euo pipefail
  git fetch origin --prune
  git branch -f main origin/main
  echo "main synced to origin/main"

git-wt-start branch name='':
  #!/usr/bin/env bash
  set -euo pipefail
  branch="{{branch}}"
  name="{{name}}"
  if [[ -z "$name" ]]; then
    name="${branch//\//-}"
  fi
  path=".worktrees/$name"
  git fetch origin --prune
  git branch -f main origin/main
  git worktree add "$path" -b "$branch" main
  echo "created worktree: $path ($branch)"

git-wt-clean path branch='':
  #!/usr/bin/env bash
  set -euo pipefail
  path="{{path}}"
  branch="{{branch}}"
  git worktree remove "$path"
  git worktree prune
  if [[ -n "$branch" ]]; then
    git branch -d "$branch"
  fi
  git fetch origin --prune
  git branch -f main origin/main
  echo "cleaned worktree: $path"
