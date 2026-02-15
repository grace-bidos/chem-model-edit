#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
usage: scripts/api/regenerate_contract.sh

Exports FastAPI OpenAPI schema and regenerates the TS API client artifacts.
Then prints focused diff hints for review.
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

echo "[1/2] Exporting OpenAPI schema"
PYTHONPATH=apps/api uv run --project apps/api python apps/api/scripts/export_openapi.py

echo "[2/2] Regenerating API client"
pnpm -C packages/api-client run generate

tracked_files=(
  "packages/api-client/openapi/openapi.json"
  "packages/api-client/src/generated/schema.ts"
)

echo
echo "Diff hint:"
echo "  git diff -- ${tracked_files[*]}"
echo

if git diff --quiet -- "${tracked_files[@]}"; then
  echo "No contract changes detected in generated artifacts."
else
  echo "Generated artifacts changed. Review and commit them in the same PR."
  git status --short -- "${tracked_files[@]}"
fi
