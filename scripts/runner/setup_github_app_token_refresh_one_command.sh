#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
One-command setup for GitHub App token refresh automation.

Usage:
  scripts/runner/setup_github_app_token_refresh_one_command.sh [options]

Options:
  --app-id <id>                  Required.
  --installation-id <id>         Required.
  --private-key-file <path>      Required.
  --token-file <path>            Optional. Default: /etc/chem-model-edit/gh_runner_token
  --status-file <path>           Optional. Default: /var/lib/chem-model-edit/github-app-token-refresh-status.json
  --interval <duration>          Optional. Default: 30min
  --service-name <name>          Optional. Default: chem-github-app-token-refresh
  --api-base <url>               Optional. Default: https://api.github.com
  --dry-run                      Print actions only.
  -h, --help                     Show help.
USAGE
}

app_id=""
installation_id=""
private_key_file=""
token_file="/etc/chem-model-edit/gh_runner_token"
status_file="/var/lib/chem-model-edit/github-app-token-refresh-status.json"
interval="30min"
service_name="chem-github-app-token-refresh"
api_base="https://api.github.com"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --app-id) app_id="${2:-}"; shift 2 ;;
    --installation-id) installation_id="${2:-}"; shift 2 ;;
    --private-key-file) private_key_file="${2:-}"; shift 2 ;;
    --token-file) token_file="${2:-}"; shift 2 ;;
    --status-file) status_file="${2:-}"; shift 2 ;;
    --interval) interval="${2:-}"; shift 2 ;;
    --service-name) service_name="${2:-}"; shift 2 ;;
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

for tool in bash sudo; do
  command -v "$tool" >/dev/null 2>&1 || { echo "missing command: $tool" >&2; exit 1; }
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
install_script="${script_dir}/install_github_app_token_refresh_timer.sh"
if [[ ! -x "$install_script" ]]; then
  echo "missing executable: $install_script" >&2
  exit 1
fi

install_args=(
  --app-id "$app_id"
  --installation-id "$installation_id"
  --private-key-file "$private_key_file"
  --token-file "$token_file"
  --status-file "$status_file"
  --interval "$interval"
  --service-name "$service_name"
  --api-base "$api_base"
)

if [[ "$dry_run" -eq 1 ]]; then
  "$install_script" "${install_args[@]}" --dry-run
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

"$install_script" "${install_args[@]}"

sudo systemctl status "${service_name}.timer" --no-pager
sudo systemctl status "${service_name}.service" --no-pager || true
