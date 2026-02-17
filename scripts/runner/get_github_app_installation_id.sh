#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Get GitHub App installation ID for a repository using App JWT.

Usage:
  scripts/runner/get_github_app_installation_id.sh \
    --app-id 123456 \
    --private-key-file /etc/chem-model-edit/github-app.pem \
    --owner grace-bidos \
    --repo chem-model-edit

Options:
  --app-id <id>                  Required.
  --private-key-file <path>      Required.
  --owner <owner>                Required.
  --repo <repo>                  Required.
  --api-base <url>               Optional. Default: https://api.github.com
  --json                          Print full JSON response.
  -h, --help                     Show help.
USAGE
}

die() {
  printf 'Error: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "missing command: $cmd"
}

require_value() {
  local flag="$1"
  local value="${2-}"
  [[ -n "$value" && "${value:0:1}" != "-" ]] || die "missing value for ${flag}"
}

b64url() {
  openssl base64 -A | tr '+/' '-_' | tr -d '='
}

app_id=""
private_key_file=""
owner=""
repo=""
api_base="https://api.github.com"
output_json=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --app-id)
      require_value "$1" "${2-}"
      app_id="$2"
      shift 2
      ;;
    --private-key-file)
      require_value "$1" "${2-}"
      private_key_file="$2"
      shift 2
      ;;
    --owner)
      require_value "$1" "${2-}"
      owner="$2"
      shift 2
      ;;
    --repo)
      require_value "$1" "${2-}"
      repo="$2"
      shift 2
      ;;
    --api-base)
      require_value "$1" "${2-}"
      api_base="$2"
      shift 2
      ;;
    --json)
      output_json=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage
      die "unknown argument: $1"
      ;;
  esac
done

[[ -n "$app_id" ]] || { usage; die "--app-id is required"; }
[[ -n "$private_key_file" ]] || { usage; die "--private-key-file is required"; }
[[ -n "$owner" ]] || { usage; die "--owner is required"; }
[[ -n "$repo" ]] || { usage; die "--repo is required"; }

for tool in openssl curl jq mktemp date; do
  require_cmd "$tool"
done

signing_key_file="$private_key_file"
tmp_key=""
cleanup() {
  [[ -n "$tmp_key" && -f "$tmp_key" ]] && rm -f "$tmp_key"
}
trap cleanup EXIT

if [[ ! -r "$private_key_file" ]]; then
  if command -v sudo >/dev/null 2>&1 && sudo -n test -r "$private_key_file" 2>/dev/null; then
    tmp_key="$(mktemp)"
    sudo cat "$private_key_file" > "$tmp_key"
    chmod 600 "$tmp_key"
    signing_key_file="$tmp_key"
  else
    die "private key is not readable: $private_key_file (run with sudo auth first: sudo -v)"
  fi
fi

now="$(date +%s)"
iat="$((now - 60))"
exp="$((now + 540))"

header='{"alg":"RS256","typ":"JWT"}'
payload="$(printf '{"iat":%s,"exp":%s,"iss":"%s"}' "$iat" "$exp" "$app_id")"
unsigned_token="$(printf '%s' "$header" | b64url).$(printf '%s' "$payload" | b64url)"
signature="$(printf '%s' "$unsigned_token" | openssl dgst -sha256 -sign "$signing_key_file" -binary | b64url)"
jwt="${unsigned_token}.${signature}"

response_and_code="$(curl -fsS -w '\n%{http_code}' \
  -H "Accept: application/vnd.github+json" \
  -H "Authorization: Bearer ${jwt}" \
  "${api_base}/repos/${owner}/${repo}/installation" || true)"

http_code="${response_and_code##*$'\n'}"
response="${response_and_code%$'\n'*}"

if [[ -z "$http_code" || ! "$http_code" =~ ^[0-9]{3}$ ]]; then
  die "failed to query ${api_base}/repos/${owner}/${repo}/installation"
fi

if [[ "$http_code" == "404" ]]; then
  printf 'Repository installation not found (404).\n' >&2
  printf 'Check that the app is installed on %s/%s.\n' "$owner" "$repo" >&2
  jq -r '.' <<<"$response" >&2 || printf '%s\n' "$response" >&2
  exit 1
fi

if [[ ! "$http_code" =~ ^2[0-9][0-9]$ ]]; then
  printf 'GitHub API error HTTP %s.\n' "$http_code" >&2
  jq -r '.' <<<"$response" >&2 || printf '%s\n' "$response" >&2
  exit 1
fi

if [[ "$output_json" -eq 1 ]]; then
  printf '%s\n' "$response"
  exit 0
fi

installation_id="$(jq -r '.id // empty' <<<"$response")"
if [[ -z "$installation_id" ]]; then
  printf 'Failed to parse installation id.\n' >&2
  jq -r '.' <<<"$response" >&2 || true
  exit 1
fi

printf '%s\n' "$installation_id"
