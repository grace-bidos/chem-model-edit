#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Verify self-hosted runner supervisor + GitHub App token refresh setup.

Usage:
  scripts/runner/verify_self_hosted_runner_setup.sh [options]

Options:
  --owner <owner>                     Optional. Default: grace-bidos
  --repo <repo>                       Optional. Default: chem-model-edit
  --token-file <path>                 Optional. Default: /etc/chem-model-edit/gh_runner_token
  --status-file <path>                Optional. Default: /var/lib/chem-model-edit/github-app-token-refresh-status.json
  --refresh-timer <name>              Optional. Default: chem-github-app-token-refresh.timer
  --refresh-service <name>            Optional. Default: chem-github-app-token-refresh.service
  --supervisor-service <name>         Optional. Default: chem-runner-pool-supervisor.service
  --no-gh-check                       Skip GitHub API runner inventory check.
  -h, --help                          Show help.
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

owner="grace-bidos"
repo="chem-model-edit"
token_file="/etc/chem-model-edit/gh_runner_token"
status_file="/var/lib/chem-model-edit/github-app-token-refresh-status.json"
refresh_timer="chem-github-app-token-refresh.timer"
refresh_service="chem-github-app-token-refresh.service"
supervisor_service="chem-runner-pool-supervisor.service"
do_gh_check=1

while [[ $# -gt 0 ]]; do
  case "$1" in
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
    --token-file)
      require_value "$1" "${2-}"
      token_file="$2"
      shift 2
      ;;
    --status-file)
      require_value "$1" "${2-}"
      status_file="$2"
      shift 2
      ;;
    --refresh-timer)
      require_value "$1" "${2-}"
      refresh_timer="$2"
      shift 2
      ;;
    --refresh-service)
      require_value "$1" "${2-}"
      refresh_service="$2"
      shift 2
      ;;
    --supervisor-service)
      require_value "$1" "${2-}"
      supervisor_service="$2"
      shift 2
      ;;
    --no-gh-check)
      do_gh_check=0
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

for tool in sudo stat systemctl journalctl dirname ls cat; do
  require_cmd "$tool"
done

if [[ "$do_gh_check" -eq 1 ]]; then
  require_cmd gh
fi

echo "==> sudo preflight"
sudo -v

echo
echo "==> Systemd status (${refresh_timer})"
sudo systemctl status "$refresh_timer" --no-pager -l

echo
echo "==> Systemd status (${refresh_service})"
sudo systemctl status "$refresh_service" --no-pager -l || true

echo
echo "==> Systemd status (${supervisor_service})"
sudo systemctl status "$supervisor_service" --no-pager -l

echo
echo "==> Recent logs (${refresh_service})"
sudo journalctl -u "$refresh_service" -n 50 --no-pager || true

echo
echo "==> Token file metadata (${token_file})"
sudo stat "$token_file"

echo
echo "==> Refresh status file (${status_file})"
if sudo test -f "$status_file"; then
  sudo cat "$status_file"
else
  echo "status file not found; listing parent directory"
  sudo ls -la "$(dirname "$status_file")"
fi

if [[ "$do_gh_check" -eq 1 ]]; then
  echo
  echo "==> GitHub runner inventory (${owner}/${repo})"
  gh api "repos/${owner}/${repo}/actions/runners" --jq '.runners[] | {name,status,busy,labels:[.labels[].name]}'
fi

echo
echo "Verification completed."
