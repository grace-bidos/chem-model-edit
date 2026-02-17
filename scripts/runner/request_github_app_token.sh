#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Request a short-lived GitHub App installation token.

Usage:
  scripts/runner/request_github_app_token.sh \
    --app-id <id> \
    --installation-id <id> \
    --private-key-file <path>

Options:
  --app-id <id>                  Required.
  --installation-id <id>         Required.
  --private-key-file <path>      Required. PEM private key for the app.
  --api-base <url>               Optional. Default: https://api.github.com
  --token-only                   Print only token (default output).
  --json                         Print raw JSON payload from GitHub API.
  -h, --help                     Show help.
EOF
}

b64url() {
  openssl base64 -A | tr '+/' '-_' | tr -d '='
}

app_id=""
installation_id=""
private_key_file=""
api_base="https://api.github.com"
output_mode="token"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --app-id) app_id="${2:-}"; shift 2 ;;
    --installation-id) installation_id="${2:-}"; shift 2 ;;
    --private-key-file) private_key_file="${2:-}"; shift 2 ;;
    --api-base) api_base="${2:-}"; shift 2 ;;
    --token-only) output_mode="token"; shift ;;
    --json) output_mode="json"; shift ;;
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

for tool in openssl curl jq; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "missing command: $tool" >&2
    exit 1
  }
done

now="$(date +%s)"
iat="$((now - 60))"
exp="$((now + 540))"

header='{"alg":"RS256","typ":"JWT"}'
payload="$(printf '{"iat":%s,"exp":%s,"iss":"%s"}' "$iat" "$exp" "$app_id")"
unsigned_token="$(printf '%s' "$header" | b64url).$(printf '%s' "$payload" | b64url)"
signature="$(printf '%s' "$unsigned_token" | openssl dgst -sha256 -sign "$private_key_file" -binary | b64url)"
jwt="${unsigned_token}.${signature}"

response="$(curl -fsSL \
  -X POST \
  -H "Accept: application/vnd.github+json" \
  -H "Authorization: Bearer ${jwt}" \
  "${api_base}/app/installations/${installation_id}/access_tokens")"

if [[ "$output_mode" == "json" ]]; then
  printf '%s\n' "$response"
  exit 0
fi

token="$(jq -r '.token // empty' <<<"$response")"
if [[ -z "$token" ]]; then
  echo "Failed to get installation token." >&2
  jq -r '.' <<<"$response" >&2 || true
  exit 1
fi
printf '%s\n' "$token"
