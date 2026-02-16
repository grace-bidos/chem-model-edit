#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  scripts/runner/cleanup_vm_temp_credentials.sh [options]

Safely remove local temporary VM credential artifacts from ignored paths.

Options:
  --path <path>   Add an extra target path (repeatable).
  --dry-run       Show what would be removed without deleting.
  --force         Delete without confirmation prompt.
  -h, --help      Show help.

Default target paths:
  .tmp/vm-validation-secrets
  .tmp/vm-jit-config
  .secrets/vm-validation

Safety guard:
  Only paths under .tmp/ or .secrets/ in this repository are allowed.
USAGE
}

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DRY_RUN=0
FORCE=0

declare -a TARGETS=(
  ".tmp/vm-validation-secrets"
  ".tmp/vm-jit-config"
  ".secrets/vm-validation"
)

declare -a NORMALIZED_TARGETS=()
declare -a REMOVABLE_PATHS=()

die() {
  echo "error: $*" >&2
  exit 1
}

normalize_path() {
  local raw="$1"
  local abs

  if [[ "$raw" = /* ]]; then
    abs="$raw"
  else
    abs="$ROOT_DIR/$raw"
  fi

  if command -v realpath >/dev/null 2>&1; then
    realpath -m "$abs"
    return
  fi

  if command -v readlink >/dev/null 2>&1; then
    readlink -m "$abs"
    return
  fi

  die "realpath or readlink is required for safe path normalization"
}

is_allowed_path() {
  local path="$1"
  case "$path" in
    "$ROOT_DIR"/.tmp/*|"$ROOT_DIR"/.secrets/*)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --path)
      [[ $# -gt 1 ]] || die "--path requires a value"
      TARGETS+=("$2")
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

for target in "${TARGETS[@]}"; do
  normalized="$(normalize_path "$target")"
  NORMALIZED_TARGETS+=("$normalized")
  if ! is_allowed_path "$normalized"; then
    die "refusing to touch path outside .tmp/ or .secrets/: $target -> $normalized"
  fi

done

for path in "${NORMALIZED_TARGETS[@]}"; do
  if [[ -e "$path" ]]; then
    REMOVABLE_PATHS+=("$path")
  fi
done

if [[ ${#REMOVABLE_PATHS[@]} -eq 0 ]]; then
  echo "No local VM credential artifacts found."
  exit 0
fi

echo "Cleanup targets:"
for path in "${REMOVABLE_PATHS[@]}"; do
  echo "- $path"
done

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry-run enabled. No files were removed."
  exit 0
fi

if [[ "$FORCE" -ne 1 ]]; then
  read -r -p "Proceed with deletion? [y/N] " response
  if [[ ! "$response" =~ ^[Yy]$ ]]; then
    echo "Aborted."
    exit 1
  fi
fi

for path in "${REMOVABLE_PATHS[@]}"; do
  rm -rf -- "$path"
done

echo "Cleanup complete."
