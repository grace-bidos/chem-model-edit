#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Validate GitHub App JWT/key pair and show app identity.

Usage:
  scripts/runner/check_github_app_jwt.sh --app-id 123456 --private-key-file /path/to/key.pem
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
api_base="https://api.github.com"

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
    --api-base)
      require_value "$1" "${2-}"
      api_base="$2"
      shift 2
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
    die "private key is not readable: $private_key_file (run: sudo -v)"
  fi
fi

now="$(date +%s)"
iat="$((now - 60))"
exp="$((now + 540))"
header='{"alg":"RS256","typ":"JWT"}'
payload="$(printf '{"iat":%s,"exp":%s,"iss":"%s"}' "$iat" "$exp" "$app_id")"
unsigned="$(printf '%s' "$header" | b64url).$(printf '%s' "$payload" | b64url)"
sig="$(printf '%s' "$unsigned" | openssl dgst -sha256 -sign "$signing_key_file" -binary | b64url)"
jwt="${unsigned}.${sig}"

resp_and_code="$(curl -fsS -w '\n%{http_code}' \
  -H "Accept: application/vnd.github+json" \
  -H "Authorization: Bearer ${jwt}" \
  "${api_base}/app" || true)"

code="${resp_and_code##*$'\n'}"
resp="${resp_and_code%$'\n'*}"

if [[ -z "$code" || ! "$code" =~ ^[0-9]{3}$ ]]; then
  die "failed to query ${api_base}/app"
fi

if [[ ! "$code" =~ ^2[0-9][0-9]$ ]]; then
  printf 'JWT validation failed: HTTP %s\n' "$code" >&2
  jq -r '.' <<<"$resp" >&2 || printf '%s\n' "$resp" >&2
  exit 1
fi

jq -r '{id,slug,name,owner:.owner.login}' <<<"$resp"
