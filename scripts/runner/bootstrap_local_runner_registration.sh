#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/runner/bootstrap_local_runner_registration.sh \
    --owner grace-bidos \
    --repo chem-model-edit \
    [--runner-home /opt/actions-runner/actions-runner] \
    [--repo-url https://github.com/grace-bidos/chem-model-edit] \
    [--runner-name <name>] \
    [--work-folder <folder>] \
    [--runner-group <name>] \
    [--labels "self-hosted,linux,x64,chem-trusted-pr"] \
    [--service-user <user>] \
    [--dry-run]

Purpose:
  Idempotently repair local runner registration when .runner is missing/invalid
  or the runner was deleted remotely. Re-registers with GitHub, then reinstalls
  and restarts the runner systemd service.

Required:
  --owner <org-or-user>
  --repo <repo-name>

Defaults:
  --runner-home /opt/actions-runner/actions-runner
  --repo-url https://github.com/<owner>/<repo>
  --runner-group Default
  --labels self-hosted,linux,x64,chem-trusted-pr
  --runner-name from existing .runner agentName, fallback: chem-base-<hostname>
  --work-folder from existing .runner workFolder, fallback: _work
  --service-user current user (or SUDO_USER when run as root)

Notes:
  - Requires: gh, jq, hostname, systemctl.
  - Preserves existing work folder; does not delete runner workspace.
  - --dry-run prints commands without mutating registration or service state.
EOF
}

log() {
  printf '[runner-bootstrap] %s\n' "$*"
}

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

require_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || die "Missing required command: $cmd"
}

quote_cmd() {
  local out=""
  local arg=""
  for arg in "$@"; do
    out+=" $(printf '%q' "$arg")"
  done
  printf '%s' "${out# }"
}

run_cmd() {
  if [[ "$dry_run" == "true" ]]; then
    printf '[dry-run] %s\n' "$(quote_cmd "$@")"
    return 0
  fi
  "$@"
}

run_privileged() {
  if [[ "$EUID" -eq 0 ]]; then
    run_cmd "$@"
    return
  fi
  run_cmd sudo "$@"
}

owner=""
repo=""
runner_home="/opt/actions-runner/actions-runner"
repo_url=""
runner_name=""
work_folder=""
runner_group="Default"
labels="self-hosted,linux,x64,chem-trusted-pr"
service_user=""
dry_run="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --owner) owner="${2:-}"; shift 2 ;;
    --repo) repo="${2:-}"; shift 2 ;;
    --runner-home) runner_home="${2:-}"; shift 2 ;;
    --repo-url) repo_url="${2:-}"; shift 2 ;;
    --runner-name) runner_name="${2:-}"; shift 2 ;;
    --work-folder) work_folder="${2:-}"; shift 2 ;;
    --runner-group) runner_group="${2:-}"; shift 2 ;;
    --labels) labels="${2:-}"; shift 2 ;;
    --service-user) service_user="${2:-}"; shift 2 ;;
    --dry-run) dry_run="true"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$owner" ]] || die "--owner is required."
[[ -n "$repo" ]] || die "--repo is required."
[[ -n "$runner_home" ]] || die "--runner-home must not be empty."
[[ -n "$runner_group" ]] || die "--runner-group must not be empty."
[[ -n "$labels" ]] || die "--labels must not be empty."

if [[ -z "$repo_url" ]]; then
  repo_url="https://github.com/${owner}/${repo}"
fi

if [[ -z "$service_user" ]]; then
  if [[ "$EUID" -eq 0 && -n "${SUDO_USER:-}" ]]; then
    service_user="${SUDO_USER}"
  else
    service_user="$(id -un)"
  fi
fi

require_cmd gh
require_cmd jq
require_cmd hostname
require_cmd systemctl

[[ -d "$runner_home" ]] || die "Runner home not found: $runner_home"
[[ -x "$runner_home/config.sh" ]] || die "Missing executable: $runner_home/config.sh"
[[ -x "$runner_home/svc.sh" ]] || die "Missing executable: $runner_home/svc.sh"

if ! gh auth status -h github.com >/dev/null 2>&1; then
  die "gh is not authenticated for github.com. Run: gh auth login"
fi

if [[ "$EUID" -ne 0 && "$dry_run" != "true" ]]; then
  require_cmd sudo
  if ! sudo -n true >/dev/null 2>&1; then
    die "sudo without prompt is required to manage systemd service. Re-run as root or configure passwordless sudo."
  fi
fi

runner_config="${runner_home}/.runner"
existing_runner_name=""
existing_work_folder=""
local_registration_state="missing"

if [[ -f "$runner_config" ]]; then
  if jq -e . "$runner_config" >/dev/null 2>&1; then
    existing_runner_name="$(jq -r '.agentName // empty' "$runner_config")"
    existing_work_folder="$(jq -r '.workFolder // empty' "$runner_config")"
    if [[ -n "$existing_runner_name" ]]; then
      local_registration_state="present"
    else
      local_registration_state="invalid"
    fi
  else
    local_registration_state="invalid"
  fi
fi

if [[ -z "$runner_name" ]]; then
  if [[ -n "$existing_runner_name" ]]; then
    runner_name="$existing_runner_name"
  else
    runner_name="chem-base-$(hostname -s)"
  fi
fi

if [[ -z "$work_folder" ]]; then
  if [[ -n "$existing_work_folder" ]]; then
    work_folder="$existing_work_folder"
  else
    work_folder="_work"
  fi
fi

log "Runner home: $runner_home"
log "Runner name: $runner_name"
log "Work folder: $work_folder"
log "Runner group: $runner_group"
log "Labels: $labels"
log "Local registration state: $local_registration_state"

remote_runner_id="$(gh api "repos/${owner}/${repo}/actions/runners?per_page=100" | jq -r --arg n "$runner_name" '.runners[]? | select(.name == $n) | .id' | head -n 1)"
needs_reregister="false"
reason=""

case "$local_registration_state" in
  missing)
    needs_reregister="true"
    reason="Local .runner config is missing."
    ;;
  invalid)
    needs_reregister="true"
    reason="Local .runner config is invalid."
    ;;
  present)
    if [[ -z "$remote_runner_id" ]]; then
      needs_reregister="true"
      reason="Runner is configured locally but not registered in GitHub."
    fi
    ;;
  *)
    die "Unexpected local state: $local_registration_state"
    ;;
esac

if [[ "$needs_reregister" == "true" ]]; then
  log "Re-registration required: $reason"

  run_privileged "$runner_home/svc.sh" stop || true
  run_privileged "$runner_home/svc.sh" uninstall || true

  run_cmd rm -f \
    "$runner_home/.runner" \
    "$runner_home/.credentials" \
    "$runner_home/.credentials_rsaparams" \
    "$runner_home/.service"

  run_cmd mkdir -p "$runner_home/$work_folder"

  reg_token="DRY_RUN_TOKEN"
  if [[ "$dry_run" != "true" ]]; then
    reg_token="$(gh api --method POST "repos/${owner}/${repo}/actions/runners/registration-token" --jq .token)"
  else
    log "[dry-run] Would request fresh registration token via gh api."
  fi
  [[ -n "$reg_token" ]] || die "Failed to obtain registration token."

  run_cmd "$runner_home/config.sh" \
    --unattended \
    --url "$repo_url" \
    --token "$reg_token" \
    --name "$runner_name" \
    --runnergroup "$runner_group" \
    --work "$work_folder" \
    --labels "$labels" \
    --replace

  run_privileged "$runner_home/svc.sh" install "$service_user"
  run_privileged "$runner_home/svc.sh" start
  log "Runner re-registration bootstrap completed."
  exit 0
fi

log "Registration is valid in both local config and GitHub."
if [[ -f "$runner_home/.service" ]]; then
  log "Restarting existing runner service."
  run_privileged "$runner_home/svc.sh" stop || true
  run_privileged "$runner_home/svc.sh" start
else
  log "Runner service metadata is missing; reinstalling service."
  run_privileged "$runner_home/svc.sh" install "$service_user"
  run_privileged "$runner_home/svc.sh" start
fi

log "Bootstrap check completed with no re-registration needed."
