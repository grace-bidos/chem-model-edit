#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SECRETS_FILE="${DEV_SECRETS_FILE:-$ROOT_DIR/.tmp/dev-secrets.env}"

usage() {
  cat <<'USAGE'
Usage: bootstrap_dev_secrets_file.sh [--file <path>]

Create or repair a local-only dev secrets file with secure permissions.
Existing values are preserved; missing keys are appended as empty entries.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --file)
      shift
      [[ $# -gt 0 ]] || { echo "error: --file requires a value" >&2; exit 1; }
      SECRETS_FILE="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

mkdir -p "$(dirname "$SECRETS_FILE")"
chmod 700 "$(dirname "$SECRETS_FILE")" || true

if [[ ! -f "$SECRETS_FILE" ]]; then
  : > "$SECRETS_FILE"
fi
chmod 600 "$SECRETS_FILE"

ensure_key() {
  local key="$1"
  if ! rg -q "^${key}=" "$SECRETS_FILE"; then
    printf '%s=\n' "$key" >> "$SECRETS_FILE"
  fi
}

# Canonical local secrets keys
ensure_key "MODAL_API_BASE_URL"
ensure_key "MODAL_PROXY_KEY"
ensure_key "MODAL_PROXY_SECRET"
ensure_key "CLERK_TEST_JWT"
ensure_key "RUNTIME_COMMAND_SUBMIT_URL"
ensure_key "RUNTIME_SERVICE_AUTH_BEARER_TOKEN"
ensure_key "REMOTE_VM_HOST"
ensure_key "REMOTE_VM_REPO_DIR"
ensure_key "REMOTE_VM_SSH_KEY"

echo "Dev secrets file ready: $SECRETS_FILE"
echo "Permissions: $(stat -c '%a' "$SECRETS_FILE")"
echo "Next: fill values, then run scripts/dev/apply_dev_runtime_env.sh"
