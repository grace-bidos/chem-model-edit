#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Scale local self-hosted runner pool for one repository.

Usage:
  scripts/runner/scale_local_runner_pool.sh \
    --repo grace-bidos/chem-model-edit \
    --baseline 1 \
    --max 4 \
    --target 4

Options:
  --repo <owner/repo>           Required.
  --target <n>                  Optional. Desired runner count (>=0).
                                If omitted, baseline is used.
  --baseline <n>                Optional floor for runner count (>=0).
                                Default: 1
  --max <n>                     Optional upper bound for runner count (>=1).
                                Default: 4
  --runner-user <user>          Runner OS user. Default: runner-user
  --base-dir <dir>              Runner root dir. Default: /opt/actions-runner
  --base-name <name>            Runner base name. Default: home-self-host
  --labels <csv>                Labels for new runners.
                                Default: self-hosted,Linux,X64,chem-trusted-pr
  --runner-version <ver>        Actions runner version. Default: 2.331.0
  --seed-dir <dir>              Seed runner install to copy from.
                                Default: /opt/actions-runner/actions-runner
  --force-remove-busy           Allow removing busy runners when scaling down.
  --dry-run                     Print operations only.
  -h, --help                    Show this help.

Naming:
  index 1 => <base-name>
  index n>=2 => <base-name>-<n>

Notes:
  - Requires: gh, jq, sudo, systemctl
  - New runners are configured as ephemeral.
  - Script uses repository registration tokens for add/remove operations.
  - Effective target is clamped to [baseline, max].
EOF
}

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "missing command: $1" >&2
    exit 1
  }
}

run_cmd() {
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[dry-run] $*"
    return 0
  fi
  "$@"
}

run_sudo() {
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[dry-run] sudo $*"
    return 0
  fi
  sudo "$@"
}

runner_name_for_index() {
  local idx="$1"
  if [[ "$idx" -eq 1 ]]; then
    printf '%s\n' "$BASE_NAME"
  else
    printf '%s-%s\n' "$BASE_NAME" "$idx"
  fi
}

runner_dir_for_name() {
  local name="$1"
  printf '%s/%s\n' "$BASE_DIR" "$name"
}

runner_exists_in_github() {
  local name="$1"
  gh api "repos/${REPO}/actions/runners" \
    --jq ".runners[] | select(.name==\"${name}\") | .id" \
    2>/dev/null | sed -n '1p'
}

runner_busy_in_github() {
  local name="$1"
  gh api "repos/${REPO}/actions/runners" \
    --jq ".runners[] | select(.name==\"${name}\") | .busy" \
    2>/dev/null | sed -n '1p'
}

registration_token() {
  if [[ "$DRY_RUN" == "true" ]]; then
    printf 'DRY_RUN_TOKEN\n'
    return 0
  fi
  gh api -X POST "repos/${REPO}/actions/runners/registration-token" --jq '.token'
}

remove_token() {
  if [[ "$DRY_RUN" == "true" ]]; then
    printf 'DRY_RUN_TOKEN\n'
    return 0
  fi
  gh api -X POST "repos/${REPO}/actions/runners/remove-token" --jq '.token'
}

configure_runner() {
  local name="$1"
  local dir
  dir="$(runner_dir_for_name "$name")"
  local token
  token="$(registration_token)"

  run_sudo mkdir -p "$dir"
  run_sudo chown -R "$RUNNER_USER:$RUNNER_USER" "$dir"

  if [[ ! -x "${dir}/config.sh" ]]; then
    if [[ -d "$SEED_DIR" && -x "${SEED_DIR}/config.sh" ]]; then
      run_sudo rsync -a --delete \
        --exclude ".runner" \
        --exclude ".credentials" \
        --exclude ".credentials_rsaparams" \
        --exclude ".service" \
        --exclude "_diag" \
        --exclude "_work" \
        --exclude "*.log" \
        "${SEED_DIR}/" "${dir}/"
    else
      local archive="actions-runner-linux-x64-${RUNNER_VERSION}.tar.gz"
      local url="https://github.com/actions/runner/releases/download/v${RUNNER_VERSION}/${archive}"
      run_sudo -u "$RUNNER_USER" bash -lc "cd '${dir}' && curl -fsSL '${url}' -o '${archive}' && tar xzf '${archive}' && rm -f '${archive}'"
    fi
  fi

  run_sudo -u "$RUNNER_USER" bash -lc "
    set -euo pipefail
    cd '${dir}'
    ./config.sh \
      --url 'https://github.com/${REPO}' \
      --token '${token}' \
      --name '${name}' \
      --labels '${LABELS}' \
      --work '_work' \
      --ephemeral \
      --unattended \
      --replace
  "

  run_sudo bash -lc "cd '${dir}' && ./svc.sh install '${RUNNER_USER}'"
  run_sudo bash -lc "cd '${dir}' && ./svc.sh start"
}

remove_runner() {
  local name="$1"
  local dir
  dir="$(runner_dir_for_name "$name")"

  local busy
  busy="$(runner_busy_in_github "$name" || true)"
  if [[ "$busy" == "true" && "$FORCE_REMOVE_BUSY" != "true" ]]; then
    echo "skip removing busy runner: ${name}" >&2
    return 0
  fi

  if [[ -d "$dir" && -x "${dir}/svc.sh" ]]; then
    run_sudo bash -lc "cd '${dir}' && ./svc.sh stop || true"
    run_sudo bash -lc "cd '${dir}' && ./svc.sh uninstall || true"
  fi

  if [[ -d "$dir" && -x "${dir}/config.sh" ]]; then
    local token
    token="$(remove_token)"
    run_sudo -u "$RUNNER_USER" bash -lc "
      set -euo pipefail
      cd '${dir}'
      ./config.sh remove --token '${token}' || true
    "
  fi

  local rid
  rid="$(runner_exists_in_github "$name" || true)"
  if [[ -n "$rid" ]]; then
    run_cmd gh api -X DELETE "repos/${REPO}/actions/runners/${rid}"
  fi

  run_sudo rm -rf "$dir"
}

REPO=""
TARGET=""
BASELINE="1"
MAX_PARALLEL="4"
RUNNER_USER="runner-user"
BASE_DIR="/opt/actions-runner"
BASE_NAME="home-self-host"
LABELS="self-hosted,Linux,X64,chem-trusted-pr"
RUNNER_VERSION="2.331.0"
SEED_DIR="/opt/actions-runner/actions-runner"
FORCE_REMOVE_BUSY="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) REPO="${2:-}"; shift 2 ;;
    --target) TARGET="${2:-}"; shift 2 ;;
    --baseline) BASELINE="${2:-}"; shift 2 ;;
    --max) MAX_PARALLEL="${2:-}"; shift 2 ;;
    --runner-user) RUNNER_USER="${2:-}"; shift 2 ;;
    --base-dir) BASE_DIR="${2:-}"; shift 2 ;;
    --base-name) BASE_NAME="${2:-}"; shift 2 ;;
    --labels) LABELS="${2:-}"; shift 2 ;;
    --runner-version) RUNNER_VERSION="${2:-}"; shift 2 ;;
    --seed-dir) SEED_DIR="${2:-}"; shift 2 ;;
    --force-remove-busy) FORCE_REMOVE_BUSY="true"; shift 1 ;;
    --dry-run) DRY_RUN="true"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$REPO" ]]; then
  echo "--repo is required." >&2
  usage
  exit 1
fi
if ! [[ "$BASELINE" =~ ^[0-9]+$ ]]; then
  echo "--baseline must be a non-negative integer." >&2
  exit 1
fi
if ! [[ "$MAX_PARALLEL" =~ ^[0-9]+$ ]]; then
  echo "--max must be a non-negative integer." >&2
  exit 1
fi
if [[ "$MAX_PARALLEL" -lt 1 ]]; then
  echo "--max must be >= 1." >&2
  exit 1
fi
if [[ "$BASELINE" -gt "$MAX_PARALLEL" ]]; then
  echo "--baseline must be <= --max." >&2
  exit 1
fi
if [[ -n "$TARGET" && ! "$TARGET" =~ ^[0-9]+$ ]]; then
  echo "--target must be a non-negative integer when provided." >&2
  exit 1
fi

requested_target="${TARGET:-$BASELINE}"
effective_target="$requested_target"
if [[ "$effective_target" -lt "$BASELINE" ]]; then
  effective_target="$BASELINE"
fi
if [[ "$effective_target" -gt "$MAX_PARALLEL" ]]; then
  effective_target="$MAX_PARALLEL"
fi

require_cmd gh
require_cmd jq
require_cmd sudo
require_cmd systemctl
require_cmd rsync

if [[ "$DRY_RUN" != "true" ]] && ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

echo "Scaling runner pool for ${REPO}"
echo "baseline=${BASELINE} max=${MAX_PARALLEL} requested_target=${requested_target} effective_target=${effective_target}"
echo "base_name=${BASE_NAME} base_dir=${BASE_DIR}"

# Build desired name set.
declare -A desired=()
for ((i=1; i<=effective_target; i++)); do
  n="$(runner_name_for_index "$i")"
  desired["$n"]=1
done

# Ensure desired runners exist.
for ((i=1; i<=effective_target; i++)); do
  n="$(runner_name_for_index "$i")"
  d="$(runner_dir_for_name "$n")"
  if [[ -d "$d" && -x "${d}/config.sh" && -f "${d}/.runner" ]]; then
    echo "exists (configured): ${n}"
    run_sudo bash -lc "cd '${d}' && ./svc.sh start || true"
    continue
  fi
  if [[ -d "$d" && -x "${d}/config.sh" ]]; then
    echo "reconfigure (missing .runner): ${n}"
  else
    echo "create: ${n}"
  fi
  configure_runner "$n"
done

# Scale down extra managed runners matching base pattern.
mapfile -t all_names < <(gh api "repos/${REPO}/actions/runners" --jq '.runners[].name' 2>/dev/null || true)
for n in "${all_names[@]}"; do
  if [[ "$n" == "$BASE_NAME" || "$n" =~ ^${BASE_NAME}-[0-9]+$ ]]; then
    if [[ -z "${desired[$n]+x}" ]]; then
      echo "remove: ${n}"
      remove_runner "$n"
    fi
  fi
done

echo "Done."
