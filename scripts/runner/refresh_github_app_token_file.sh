#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Refresh GitHub App installation token and rotate token file atomically.

Usage:
  scripts/runner/refresh_github_app_token_file.sh \
    --app-id <id> \
    --installation-id <id> \
    --private-key-file <path> \
    --token-file /etc/chem-model-edit/gh_runner_token

Options:
  --app-id <id>                  Required.
  --installation-id <id>         Required.
  --private-key-file <path>      Required.
  --token-file <path>            Optional. Default: /etc/chem-model-edit/gh_runner_token
  --status-file <path>           Optional. Default: /var/lib/chem-model-edit/github-app-token-refresh-status.json
  --api-base <url>               Optional. Default: https://api.github.com
  --dry-run                      Print actions only.
  -h, --help                     Show help.
USAGE
}

as_root() {
  if [[ "$(id -u)" -eq 0 ]]; then
    "$@"
  else
    sudo "$@"
  fi
}

app_id=""
installation_id=""
private_key_file=""
token_file="/etc/chem-model-edit/gh_runner_token"
status_file="/var/lib/chem-model-edit/github-app-token-refresh-status.json"
api_base="https://api.github.com"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --app-id) app_id="${2:-}"; shift 2 ;;
    --installation-id) installation_id="${2:-}"; shift 2 ;;
    --private-key-file) private_key_file="${2:-}"; shift 2 ;;
    --token-file) token_file="${2:-}"; shift 2 ;;
    --status-file) status_file="${2:-}"; shift 2 ;;
    --api-base) api_base="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$app_id" || -z "$installation_id" || -z "$private_key_file" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if [[ ! -f "$private_key_file" ]]; then
  echo "Private key file not found: $private_key_file" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
request_script="${script_dir}/request_github_app_token.sh"
if [[ ! -x "$request_script" ]]; then
  echo "missing executable: $request_script" >&2
  exit 1
fi

for tool in jq mktemp date; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "missing command: $tool" >&2
    exit 1
  }
done

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] request installation token (app_id=${app_id}, installation_id=${installation_id})"
  echo "[dry-run] rotate token file atomically: ${token_file}"
  echo "[dry-run] update status file: ${status_file}"
  exit 0
fi

if [[ "$(id -u)" -ne 0 ]] && ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

response_json="$($request_script \
  --app-id "$app_id" \
  --installation-id "$installation_id" \
  --private-key-file "$private_key_file" \
  --api-base "$api_base" \
  --json)"

token="$(jq -r '.token // empty' <<<"$response_json")"
expires_at="$(jq -r '.expires_at // empty' <<<"$response_json")"
if [[ -z "$token" || -z "$expires_at" ]]; then
  echo "Failed to parse token response." >&2
  jq -r '.' <<<"$response_json" >&2 || true
  exit 1
fi

token_dir="$(dirname "$token_file")"
status_dir="$(dirname "$status_file")"

local_tmp="$(mktemp)"
trap 'rm -f "$local_tmp"' EXIT
printf '%s\n' "$token" > "$local_tmp"

as_root install -d -m 0700 "$token_dir"
as_root install -m 0600 "$local_tmp" "${token_file}.new"
as_root mv -f "${token_file}.new" "$token_file"

now_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
status_tmp="$(mktemp)"
trap 'rm -f "$local_tmp" "$status_tmp"' EXIT
jq -n \
  --arg refreshed_at "$now_utc" \
  --arg expires_at "$expires_at" \
  --arg token_file "$token_file" \
  '{refreshed_at:$refreshed_at,expires_at:$expires_at,token_file:$token_file}' > "$status_tmp"

as_root install -d -m 0755 "$status_dir"
as_root install -m 0644 "$status_tmp" "${status_file}.new"
as_root mv -f "${status_file}.new" "$status_file"

echo "Token refreshed: ${token_file} (expires_at=${expires_at})"
