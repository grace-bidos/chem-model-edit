#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  RUNNER_OWNER=<org-or-user> \
  RUNNER_REPO=<repo-name> \
  RUNNER_LABELS="self-hosted,linux,x64,chem-trusted-pr" \
  RUNNER_GROUP="Default" \
  scripts/runner/recover_base_runner.sh [--dry-run] [--runner-home /opt/actions-runner/actions-runner]

Required environment variables:
  RUNNER_OWNER   GitHub org or user owner.
  RUNNER_REPO    GitHub repository name.
  RUNNER_LABELS  Comma-separated labels for runner registration.
  RUNNER_GROUP   Runner group name.

Optional environment variables:
  RUNNER_NAME           Runner name (default: local-hostname-timestamp).
  RUNNER_WORKDIR        Runner work folder (default: _work).
  RUNNER_SERVICE_USER   User for `svc.sh install` (default: current user).

Options:
  --runner-home <path>  Base runner install directory.
                        Default: /opt/actions-runner/actions-runner
  --dry-run             Print the actions without executing.
  -h, --help            Show this help.

What this does:
  1) Acquire remove and registration tokens from GitHub.
  2) Stop/uninstall existing service if present.
  3) Remove old runner config (best effort).
  4) Reconfigure runner with required labels/group.
  5) Install/start service.
USAGE
}

runner_home="/opt/actions-runner/actions-runner"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --runner-home) runner_home="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

required_env=(RUNNER_OWNER RUNNER_REPO RUNNER_LABELS RUNNER_GROUP)
for key in "${required_env[@]}"; do
  if [[ -z "${!key:-}" ]]; then
    echo "Missing required env var: ${key}" >&2
    usage
    exit 1
  fi
done

runner_name="${RUNNER_NAME:-$(hostname)-$(date +%Y%m%d-%H%M%S)}"
runner_workdir="${RUNNER_WORKDIR:-_work}"
runner_service_user="${RUNNER_SERVICE_USER:-$(id -un)}"
repo_slug="${RUNNER_OWNER}/${RUNNER_REPO}"

for tool in gh jq; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing required tool: $tool" >&2
    exit 2
  fi
done

run_privileged() {
  if [[ "$dry_run" -eq 1 ]]; then
    echo "[dry-run] sudo $*"
    return 0
  fi
  sudo "$@"
}

if [[ ! -d "$runner_home" ]]; then
  echo "Runner home not found: $runner_home" >&2
  exit 1
fi

if [[ ! -x "$runner_home/config.sh" || ! -x "$runner_home/svc.sh" ]]; then
  echo "Expected config.sh and svc.sh under: $runner_home" >&2
  exit 1
fi

run_cmd() {
  if [[ "$dry_run" -eq 1 ]]; then
    echo "[dry-run] $*"
  else
    "$@"
  fi
}

get_token() {
  local endpoint="$1"
  if [[ "$dry_run" -eq 1 ]]; then
    echo "DRY_RUN_TOKEN"
    return 0
  fi

  gh api --method POST "repos/${repo_slug}/actions/runners/${endpoint}" | jq -r '.token'
}

echo "Recovering base runner for ${repo_slug}"
echo "Runner home: ${runner_home}"
echo "Runner name: ${runner_name}"
echo "Runner group: ${RUNNER_GROUP}"
echo "Runner labels: ${RUNNER_LABELS}"
if [[ "$dry_run" -eq 1 ]]; then
  echo "Dry-run mode is ON (no changes executed)."
fi

remove_token="$(get_token remove-token)"
registration_token="$(get_token registration-token)"

if [[ -z "$remove_token" || "$remove_token" == "null" ]]; then
  echo "Failed to acquire remove token." >&2
  exit 1
fi
if [[ -z "$registration_token" || "$registration_token" == "null" ]]; then
  echo "Failed to acquire registration token." >&2
  exit 1
fi

cd "$runner_home"

# Refresh systemd unit cache if unit files were edited.
if command -v systemctl >/dev/null 2>&1; then
  run_privileged systemctl daemon-reload || true
fi

run_privileged ./svc.sh stop || true
run_privileged ./svc.sh uninstall || true

# Best effort removal because local config may already be partially broken.
run_cmd ./config.sh remove --unattended --token "$remove_token" || true

run_cmd ./config.sh \
  --url "https://github.com/${repo_slug}" \
  --token "$registration_token" \
  --name "$runner_name" \
  --runnergroup "$RUNNER_GROUP" \
  --labels "$RUNNER_LABELS" \
  --work "$runner_workdir" \
  --unattended \
  --replace

run_privileged ./svc.sh install "$runner_service_user"
run_privileged ./svc.sh start

echo "Recovery flow complete."
