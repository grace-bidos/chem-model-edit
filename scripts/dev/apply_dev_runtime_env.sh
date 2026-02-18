#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SECRETS_FILE="${DEV_SECRETS_FILE:-$ROOT_DIR/.tmp/dev-secrets.env}"
API_ENV_FILE="$ROOT_DIR/apps/api/.env"

usage() {
  cat <<'USAGE'
Usage: apply_dev_runtime_env.sh [--secrets-file <path>] [--api-env <path>]

Apply runtime-related local dev secrets into apps/api/.env.
This script updates only:
- RUNTIME_COMMAND_SUBMIT_URL
- RUNTIME_SERVICE_AUTH_BEARER_TOKEN
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --secrets-file)
      shift
      [[ $# -gt 0 ]] || { echo "error: --secrets-file requires a value" >&2; exit 1; }
      SECRETS_FILE="$1"
      ;;
    --api-env)
      shift
      [[ $# -gt 0 ]] || { echo "error: --api-env requires a value" >&2; exit 1; }
      API_ENV_FILE="$1"
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

[[ -f "$SECRETS_FILE" ]] || { echo "error: secrets file not found: $SECRETS_FILE" >&2; exit 1; }

set -a
# shellcheck disable=SC1090
source "$SECRETS_FILE"
set +a

: "${RUNTIME_COMMAND_SUBMIT_URL:?RUNTIME_COMMAND_SUBMIT_URL is required in $SECRETS_FILE}"
: "${RUNTIME_SERVICE_AUTH_BEARER_TOKEN:?RUNTIME_SERVICE_AUTH_BEARER_TOKEN is required in $SECRETS_FILE}"

mkdir -p "$(dirname "$API_ENV_FILE")"
[[ -f "$API_ENV_FILE" ]] || : > "$API_ENV_FILE"

upsert() {
  local key="$1"
  local value="$2"
  if rg -q "^${key}=" "$API_ENV_FILE"; then
    sed -i "s#^${key}=.*#${key}=${value}#" "$API_ENV_FILE"
  else
    printf '%s=%s\n' "$key" "$value" >> "$API_ENV_FILE"
  fi
}

upsert "RUNTIME_COMMAND_SUBMIT_URL" "$RUNTIME_COMMAND_SUBMIT_URL"
upsert "RUNTIME_SERVICE_AUTH_BEARER_TOKEN" "$RUNTIME_SERVICE_AUTH_BEARER_TOKEN"

echo "Applied runtime settings into: $API_ENV_FILE"
echo "- RUNTIME_COMMAND_SUBMIT_URL: set"
echo "- RUNTIME_SERVICE_AUTH_BEARER_TOKEN: set"
