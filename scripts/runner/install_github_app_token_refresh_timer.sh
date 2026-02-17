#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Install systemd timer/service to refresh GitHub App installation token.

Usage:
  scripts/runner/install_github_app_token_refresh_timer.sh \
    --app-id <id> \
    --installation-id <id> \
    --private-key-file <path>

Options:
  --app-id <id>                  Required.
  --installation-id <id>         Required.
  --private-key-file <path>      Required.
  --token-file <path>            Optional. Default: /etc/chem-model-edit/gh_runner_token
  --status-file <path>           Optional. Default: /var/lib/chem-model-edit/github-app-token-refresh-status.json
  --interval <systemd duration>  Optional. Default: 30min
  --service-name <name>          Optional. Default: chem-github-app-token-refresh
  --repo-root <path>             Optional. Default: current git root
  --api-base <url>               Optional. Default: https://api.github.com
  --dry-run                      Print files/commands only.
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
repo_root=""
api_base="https://api.github.com"
dry_run=0
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --app-id) app_id="${2:-}"; shift 2 ;;
    --installation-id) installation_id="${2:-}"; shift 2 ;;
    --private-key-file) private_key_file="${2:-}"; shift 2 ;;
    --token-file) token_file="${2:-}"; shift 2 ;;
    --status-file) status_file="${2:-}"; shift 2 ;;
    --interval) interval="${2:-}"; shift 2 ;;
    --service-name) service_name="${2:-}"; shift 2 ;;
    --repo-root) repo_root="${2:-}"; shift 2 ;;
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

if [[ -z "$repo_root" ]]; then
  repo_root="$(git -C "$script_dir" rev-parse --show-toplevel)"
fi

refresh_script="${repo_root}/scripts/runner/refresh_github_app_token_file.sh"
if [[ ! -x "$refresh_script" ]]; then
  echo "missing executable: $refresh_script" >&2
  exit 1
fi

service_file="/etc/systemd/system/${service_name}.service"
timer_file="/etc/systemd/system/${service_name}.timer"

service_content="[Unit]
Description=Refresh GitHub App installation token for local runner operations
After=network-online.target
Wants=network-online.target

[Service]
Type=oneshot
ExecStart=\"${refresh_script}\" --app-id \"${app_id}\" --installation-id \"${installation_id}\" --private-key-file \"${private_key_file}\" --token-file \"${token_file}\" --status-file \"${status_file}\" --api-base \"${api_base}\"
"

timer_content="[Unit]
Description=Periodic GitHub App installation token refresh

[Timer]
OnBootSec=30s
OnUnitActiveSec=${interval}
Persistent=true

[Install]
WantedBy=timers.target
"

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] write ${service_file}:"
  printf '%s\n' "$service_content"
  echo "[dry-run] write ${timer_file}:"
  printf '%s\n' "$timer_content"
  echo "[dry-run] sudo systemctl daemon-reload"
  echo "[dry-run] sudo systemctl enable --now ${service_name}.timer"
  echo "[dry-run] sudo systemctl start ${service_name}.service"
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

tmpd="$(mktemp -d)"
trap 'rm -rf "$tmpd"' EXIT

printf '%s\n' "$service_content" > "${tmpd}/service"
printf '%s\n' "$timer_content" > "${tmpd}/timer"

sudo install -m 0644 "${tmpd}/service" "$service_file"
sudo install -m 0644 "${tmpd}/timer" "$timer_file"

sudo systemctl daemon-reload
sudo systemctl enable --now "${service_name}.timer"
sudo systemctl start "${service_name}.service"
sudo systemctl status "${service_name}.timer" --no-pager
