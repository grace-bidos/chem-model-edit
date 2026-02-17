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
  uv run mypy --config-file pyproject.toml app services main.py
  popd >/dev/null

api-pyright:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run pyright --project pyrightconfig.json
  popd >/dev/null

api-pyright-diff:
  #!/usr/bin/env bash
  set -euo pipefail
  base_ref="${BASE_REF:-origin/main}"
  uv run --project apps/api python scripts/api/pyright_touched_gate.py --base-ref "$base_ref"

api-typecheck-strict:
  #!/usr/bin/env bash
  set -euo pipefail
  just api-mypy
  just api-pyright

api-security:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run bandit -q -r app services main.py
  uv run pip-audit --progress-spinner off
  popd >/dev/null

api-security-bandit:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run bandit -q -r app services main.py
  popd >/dev/null

api-security-audit:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run pip-audit --progress-spinner off
  popd >/dev/null

api-deadcode:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  uv run deptry . --extend-exclude ".*/mutants/"
  uv run vulture app services --min-confidence 90 --ignore-names cls
  popd >/dev/null

api-schemathesis-smoke:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  SCHEMATHESIS_MODE=smoke uv run bash scripts/schemathesis.sh
  popd >/dev/null

api-schemathesis-broad:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  SCHEMATHESIS_MODE=broad uv run bash scripts/schemathesis.sh
  popd >/dev/null

api-mutation-smoke:
  #!/usr/bin/env bash
  set -euo pipefail
  pushd apps/api >/dev/null
  bash scripts/mutmut_smoke.sh
  popd >/dev/null

api-quality-fast:
  #!/usr/bin/env bash
  set -euo pipefail
  just api-ruff
  just api-typecheck-strict
  just api-test

api-quality-security-hygiene:
  #!/usr/bin/env bash
  set -euo pipefail
  just api-security
  just api-deadcode

api-quality-phase1:
  #!/usr/bin/env bash
  set -euo pipefail
  just api-quality-fast
  just api-quality-security-hygiene
  just api-schemathesis-smoke
  just api-mutation-smoke

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
  uv run pytest --cov=app --cov=services --cov=main --cov-report=term-missing --cov-report=html:htmlcov --cov-report=xml:coverage.xml
  popd >/dev/null
  echo "generated: apps/api/htmlcov/index.html"
  echo "generated: apps/api/coverage.xml"

api-cov-diff:
  #!/usr/bin/env bash
  set -euo pipefail
  base_ref="${BASE_REF:-origin/main}"
  threshold="${COVERAGE_DIFF_THRESHOLD:-90}"
  if [[ ! -f apps/api/coverage.xml ]]; then
    just api-cov
  fi
  uv run --project apps/api python - "$base_ref" "$threshold" <<'PY'
  import re
  import subprocess
  import sys
  import xml.etree.ElementTree as ET
  from pathlib import Path

  base_ref = sys.argv[1]
  threshold = float(sys.argv[2])

  coverage_path = Path("apps/api/coverage.xml")
  if not coverage_path.exists():
      print(f"coverage report not found: {coverage_path}", file=sys.stderr)
      sys.exit(2)

  tree = ET.parse(coverage_path)
  root = tree.getroot()

  measurable: dict[str, set[int]] = {}
  covered: dict[str, set[int]] = {}

  for klass in root.findall(".//class"):
      filename = klass.attrib.get("filename")
      if not filename:
          continue
      key = Path(filename).as_posix()
      measurable.setdefault(key, set())
      covered.setdefault(key, set())
      for line in klass.findall("./lines/line"):
          number_text = line.attrib.get("number")
          hits_text = line.attrib.get("hits", "0")
          if number_text is None:
              continue
          line_no = int(number_text)
          hits = int(hits_text)
          measurable[key].add(line_no)
          if hits > 0:
              covered[key].add(line_no)

  diff = subprocess.check_output(
      [
          "git",
          "diff",
          "--unified=0",
          "--no-color",
          f"{base_ref}...HEAD",
          "--",
          "apps/api",
      ],
      text=True,
  )

  changed: dict[str, set[int]] = {}
  current_file: str | None = None
  new_line: int | None = None

  for raw in diff.splitlines():
      if raw.startswith("+++ b/"):
          path = raw[6:]
          current_file = path if path.endswith(".py") else None
          if current_file is not None and current_file.startswith("apps/api/"):
              current_file = current_file[len("apps/api/") :]
          continue
      if raw.startswith("@@"):
          match = re.search(r"\+(\d+)(?:,(\d+))?", raw)
          if not match:
              new_line = None
              continue
          new_line = int(match.group(1))
          continue
      if current_file is None or new_line is None:
          continue
      if raw.startswith("+") and not raw.startswith("+++"):
          changed.setdefault(current_file, set()).add(new_line)
          new_line += 1
      elif raw.startswith("-") and not raw.startswith("---"):
          continue
      elif raw.startswith(" "):
          new_line += 1

  measured_total = 0
  covered_total = 0
  missing: list[str] = []

  for filename, changed_lines in sorted(changed.items()):
      measurable_lines = measurable.get(filename, set())
      covered_lines = covered.get(filename, set())
      for line_no in sorted(changed_lines):
          if line_no not in measurable_lines:
              continue
          measured_total += 1
          if line_no in covered_lines:
              covered_total += 1
          else:
              missing.append(f"{filename}:{line_no}")

  if measured_total == 0:
      print(f"Changed-lines coverage: no measurable Python lines changed vs {base_ref}; gate passes.")
      sys.exit(0)

  percent = (covered_total / measured_total) * 100.0
  print(
      f"Changed-lines coverage: {covered_total}/{measured_total} = {percent:.2f}% "
      f"(threshold {threshold:.2f}%, base {base_ref})"
  )
  if missing:
      print("Uncovered changed lines:")
      for item in missing:
          print(f"  - {item}")

  if percent + 1e-9 < threshold:
      sys.exit(1)
  PY

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

web-test-coverage:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:coverage

web-test-a11y:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:a11y

web-test-fastcheck:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:fastcheck

web-test-mutation:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:mutation

web-test-e2e:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:e2e

web-test-e2e-install:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm -C apps/web test:e2e:install

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
  if lsof -ti tcp:6006 -sTCP:LISTEN >/dev/null 2>&1; then
    echo "port 6006 is already in use; stop the existing process before running story-catalog-check." >&2
    exit 1
  fi
  cleanup() {
    if [[ -n "${storybook_pid:-}" ]]; then
      kill "$storybook_pid" 2>/dev/null || true
      wait "$storybook_pid" 2>/dev/null || true
    fi
  }
  trap cleanup EXIT
  pnpm -C apps/web storybook --port 6006 >/tmp/storybook.log 2>&1 &
  storybook_pid=$!
  for i in {1..60}; do
    if ! kill -0 "$storybook_pid" 2>/dev/null; then
      echo "Storybook process exited before readiness check." >&2
      echo "--- storybook log ---" >&2
      cat /tmp/storybook.log >&2
      exit 1
    fi
    if curl -sf "$url/iframe.html" | grep -qi "storybook"; then
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
  node -e "const fs=require('fs'); const p='docs/graphs/web-dependency-graph.json'; const data=JSON.parse(fs.readFileSync(p,'utf8')); if (data?.summary?.optionsUsed?.baseDir) data.summary.optionsUsed.baseDir='.'; fs.writeFileSync(p, JSON.stringify(data, null, 2) + '\n');"
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

quality-quick:
  #!/usr/bin/env bash
  set -euo pipefail
  just web-lint
  just web-typecheck
  just web-test

quality-standard:
  #!/usr/bin/env bash
  set -euo pipefail
  just quality-quick
  just web-test-a11y
  pnpm exec nx run web:knip
  just web-test-fastcheck
  just web-test-coverage

quality-deep:
  #!/usr/bin/env bash
  set -euo pipefail
  just quality-standard
  just graph-web-deps
  just storybook-build
  if [[ -n "${CHROMATIC_PROJECT_TOKEN:-}" ]]; then
    just chromatic
  else
    echo "skip chromatic: CHROMATIC_PROJECT_TOKEN is not set"
  fi
  just web-test-mutation
  just web-test-e2e-install
  just web-test-e2e

typecheck:
  #!/usr/bin/env bash
  set -euo pipefail
  just web-typecheck
  just api-typecheck-strict

ci:
  #!/usr/bin/env bash
  set -euo pipefail
  pnpm exec nx run-many -t lint,typecheck,test,knip

git-main-sync:
  #!/usr/bin/env bash
  set -euo pipefail
  git fetch origin --prune
  current_branch="$(git rev-parse --abbrev-ref HEAD)"
  if [[ "$current_branch" != "main" ]]; then
    echo "git-main-sync must run from main branch (current: $current_branch)" >&2
    echo "hint: git switch main" >&2
    exit 1
  fi
  git merge --ff-only origin/main
  echo "main fast-forwarded to origin/main"

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
  git worktree add "$path" -b "$branch" origin/main
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
  echo "cleaned worktree: $path"

git-main-realign:
  #!/usr/bin/env bash
  set -euo pipefail
  git fetch origin --prune
  if git worktree list | grep -qE '\[main\]$'; then
    echo "main is currently checked out in a worktree; cannot force realign safely." >&2
    echo "hint: remove/switch that worktree, then rerun git-main-realign" >&2
    exit 1
  fi
  git branch -f main origin/main
  echo "main branch pointer realigned to origin/main"
