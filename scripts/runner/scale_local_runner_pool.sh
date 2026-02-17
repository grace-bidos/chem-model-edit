#!/usr/bin/env bash
set -euo pipefail

SCRIPT_VERSION="2"

usage() {
  cat <<'USAGE'
Scale local self-hosted runner pool for one repository.

Usage:
  scripts/runner/scale_local_runner_pool.sh \
    --repo grace-bidos/chem-model-edit \
    --min 1 \
    --max 4 \
    --target 4

Options:
  --repo <owner/repo>           Required.
  --target <n>                  Optional. Desired runner count (>=0).
                                If omitted, min is used.
  --min <n>                     Optional floor for runner count (>=0).
                                Default: 1
  --baseline <n>                Deprecated alias for --min.
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
  --inventory-out <path>        Optional machine-readable inventory JSON output.
  --status-out <path>           Optional machine-readable run status JSON output.
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
  - Effective target is clamped to [min, max].
  - Managed runner removal order is deterministic (highest index first).
USAGE
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
  if [[ "${EUID}" -eq 0 ]]; then
    if [[ "${1:-}" == "-u" ]]; then
      local target_user="$2"
      shift 2
      runuser -u "$target_user" -- "$@"
      return 0
    fi
    "$@"
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

managed_index_for_name() {
  local name="$1"
  if [[ "$name" == "$BASE_NAME" ]]; then
    printf '1\n'
    return 0
  fi

  local suffix
  suffix="${name#"$BASE_NAME"-}"
  if [[ "$suffix" == "$name" ]]; then
    return 1
  fi
  if [[ ! "$suffix" =~ ^[0-9]+$ ]]; then
    return 1
  fi
  if [[ "$suffix" -lt 2 ]]; then
    return 1
  fi
  printf '%s\n' "$suffix"
}

add_operation() {
  local action="$1"
  local runner="$2"
  local detail="${3:-}"
  OPERATION_LINES+=("$(jq -cn --arg action "$action" --arg runner "$runner" --arg detail "$detail" '{action:$action, runner:$runner, detail:$detail}')")
}

add_warning() {
  local code="$1"
  local runner="$2"
  local message="$3"
  WARNING_LINES+=("$(jq -cn --arg code "$code" --arg runner "$runner" --arg message "$message" '{code:$code, runner:$runner, message:$message}')")
}

json_array_from_jsonl() {
  local -n arr_ref=$1
  if [[ "${#arr_ref[@]}" -eq 0 ]]; then
    printf '[]\n'
    return 0
  fi
  printf '%s\n' "${arr_ref[@]}" | jq -s '.'
}

json_array_from_lines() {
  if [[ $# -eq 0 ]]; then
    printf '[]\n'
    return 0
  fi
  printf '%s\n' "$@" | jq -Rsc 'split("\n") | map(select(length > 0))'
}

write_output_file() {
  local path="$1"
  local content="$2"
  mkdir -p "$(dirname "$path")"
  printf '%s\n' "$content" >"$path"
}

github_runners_json() {
  gh api "repos/${REPO}/actions/runners" 2>/dev/null
}

managed_name_rows_from_lines() {
  while IFS= read -r name; do
    [[ -n "$name" ]] || continue
    local idx
    if idx="$(managed_index_for_name "$name" 2>/dev/null)"; then
      printf '%s\t%s\n' "$idx" "$name"
    fi
  done
}

managed_names_from_github_sorted() {
  local order="$1"
  local sort_cmd=(sort -n -k1,1 -k2,2)
  if [[ "$order" == "desc" ]]; then
    sort_cmd=(sort -nr -k1,1 -k2,2)
  fi

  github_runners_json \
    | jq -r '.runners[]?.name // empty' \
    | managed_name_rows_from_lines \
    | "${sort_cmd[@]}" \
    | cut -f2- \
    | awk '!seen[$0]++'
}

managed_names_from_local_sorted() {
  if [[ ! -d "$BASE_DIR" ]]; then
    return 0
  fi

  find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d -printf '%f\n' \
    | managed_name_rows_from_lines \
    | sort -n -k1,1 -k2,2 \
    | cut -f2- \
    | awk '!seen[$0]++'
}

runner_exists_in_github() {
  local name="$1"
  github_runners_json | jq -r --arg name "$name" '.runners[]? | select(.name == $name) | .id' | sed -n '1p'
}

runner_busy_in_github() {
  local name="$1"
  github_runners_json | jq -r --arg name "$name" '.runners[]? | select(.name == $name) | .busy' | sed -n '1p'
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
  local mode="$2"
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

  if [[ -f "${dir}/.service" && -x "${dir}/svc.sh" ]]; then
    run_sudo bash -lc "cd '${dir}' && ./svc.sh stop || true"
    run_sudo bash -lc "cd '${dir}' && ./svc.sh uninstall || true"
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
  add_operation "$mode" "$name" "configured_and_started"
}

remove_runner() {
  local name="$1"
  local dir
  dir="$(runner_dir_for_name "$name")"

  local busy
  busy="$(runner_busy_in_github "$name" || true)"
  if [[ "$busy" == "true" && "$FORCE_REMOVE_BUSY" != "true" ]]; then
    echo "skip removing busy runner: ${name}" >&2
    add_warning "busy_runner_not_removed" "$name" "runner is busy and --force-remove-busy was not set"
    add_operation "skip_remove_busy" "$name" "busy=true"
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
  add_operation "remove" "$name" "removed_local_and_github_if_present"
}

emit_inventory() {
  local desired_json github_json local_json runners_json
  desired_json="$(json_array_from_lines "${DESIRED_NAMES[@]}")"
  github_json="$(json_array_from_lines "${GITHUB_MANAGED_NAMES[@]}")"
  local_json="$(json_array_from_lines "${LOCAL_MANAGED_NAMES[@]}")"

  declare -A union_map=()
  local name
  for name in "${DESIRED_NAMES[@]}"; do union_map["$name"]=1; done
  for name in "${GITHUB_MANAGED_NAMES[@]}"; do union_map["$name"]=1; done
  for name in "${LOCAL_MANAGED_NAMES[@]}"; do union_map["$name"]=1; done

  declare -a union_rows=()
  for name in "${!union_map[@]}"; do
    local idx
    idx="$(managed_index_for_name "$name")"
    union_rows+=("${idx}"$'\t'"${name}")
  done

  declare -a union_sorted=()
  if [[ "${#union_rows[@]}" -gt 0 ]]; then
    mapfile -t union_sorted < <(printf '%s\n' "${union_rows[@]}" | sort -n -k1,1 -k2,2 | cut -f2-)
  fi

  declare -a runner_objs=()
  for name in "${union_sorted[@]}"; do
    local idx dir
    idx="$(managed_index_for_name "$name")"
    dir="$(runner_dir_for_name "$name")"

    local desired_present=false
    local github_present=false
    local local_present=false
    local dir_exists=false
    local config_present=false
    local runner_marker_present=false

    [[ -n "${DESIRED_MAP[$name]+x}" ]] && desired_present=true
    [[ -n "${GITHUB_MAP[$name]+x}" ]] && github_present=true
    [[ -n "${LOCAL_MAP[$name]+x}" ]] && local_present=true

    if [[ -d "$dir" ]]; then
      dir_exists=true
    fi
    if [[ -x "${dir}/config.sh" ]]; then
      config_present=true
    fi
    if [[ -f "${dir}/.runner" ]]; then
      runner_marker_present=true
    fi

    runner_objs+=("$(jq -cn \
      --arg name "$name" \
      --arg dir "$dir" \
      --argjson index "$idx" \
      --argjson desired "$desired_present" \
      --argjson in_github "$github_present" \
      --argjson in_local_dirs "$local_present" \
      --argjson dir_exists "$dir_exists" \
      --argjson config_present "$config_present" \
      --argjson runner_marker_present "$runner_marker_present" \
      '{name:$name,index:$index,dir:$dir,desired:$desired,in_github:$in_github,in_local_dirs:$in_local_dirs,dir_exists:$dir_exists,config_present:$config_present,runner_marker_present:$runner_marker_present}')")
  done

  runners_json="$(json_array_from_jsonl runner_objs)"

  local content
  content="$(jq -cn \
    --arg script "scale_local_runner_pool.sh" \
    --arg script_version "$SCRIPT_VERSION" \
    --arg repo "$REPO" \
    --arg base_name "$BASE_NAME" \
    --arg base_dir "$BASE_DIR" \
    --argjson min "$MIN_POOL" \
    --argjson max "$MAX_PARALLEL" \
    --argjson requested_target "$requested_target" \
    --argjson effective_target "$effective_target" \
    --argjson dry_run "$([[ "$DRY_RUN" == "true" ]] && echo true || echo false)" \
    --argjson desired_names "$desired_json" \
    --argjson github_managed_names "$github_json" \
    --argjson local_managed_names "$local_json" \
    --argjson runners "$runners_json" \
    '{script:$script,script_version:$script_version,repo:$repo,base_name:$base_name,base_dir:$base_dir,control:{min:$min,max:$max,target_requested:$requested_target,target_effective:$effective_target},consistency:{github_runner_inventory:"eventual"},dry_run:$dry_run,desired_names:$desired_names,github_managed_names:$github_managed_names,local_managed_names:$local_managed_names,runners:$runners}')"

  write_output_file "$INVENTORY_OUT" "$content"
}

emit_status() {
  local operations_json warnings_json
  operations_json="$(json_array_from_jsonl OPERATION_LINES)"
  warnings_json="$(json_array_from_jsonl WARNING_LINES)"

  local content
  content="$(jq -cn \
    --arg script "scale_local_runner_pool.sh" \
    --arg script_version "$SCRIPT_VERSION" \
    --arg result "ok" \
    --arg repo "$REPO" \
    --arg base_name "$BASE_NAME" \
    --argjson min "$MIN_POOL" \
    --argjson max "$MAX_PARALLEL" \
    --argjson requested_target "$requested_target" \
    --argjson effective_target "$effective_target" \
    --argjson dry_run "$([[ "$DRY_RUN" == "true" ]] && echo true || echo false)" \
    --argjson desired_count "${#DESIRED_NAMES[@]}" \
    --argjson github_managed_count "${#GITHUB_MANAGED_NAMES[@]}" \
    --argjson local_managed_count "${#LOCAL_MANAGED_NAMES[@]}" \
    --argjson operations "$operations_json" \
    --argjson warnings "$warnings_json" \
    '{script:$script,script_version:$script_version,result:$result,repo:$repo,base_name:$base_name,control:{min:$min,max:$max,target_requested:$requested_target,target_effective:$effective_target},dry_run:$dry_run,counts:{desired:$desired_count,github_managed:$github_managed_count,local_managed:$local_managed_count,operations:($operations|length),warnings:($warnings|length)},operations:$operations,warnings:$warnings}')"

  write_output_file "$STATUS_OUT" "$content"
}

REPO=""
TARGET=""
MIN_POOL="1"
MIN_ALIAS=""
MAX_PARALLEL="4"
RUNNER_USER="runner-user"
BASE_DIR="/opt/actions-runner"
BASE_NAME="home-self-host"
LABELS="self-hosted,Linux,X64,chem-trusted-pr"
RUNNER_VERSION="2.331.0"
SEED_DIR="/opt/actions-runner/actions-runner"
INVENTORY_OUT=""
STATUS_OUT=""
FORCE_REMOVE_BUSY="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) REPO="${2:-}"; shift 2 ;;
    --target) TARGET="${2:-}"; shift 2 ;;
    --min) MIN_POOL="${2:-}"; shift 2 ;;
    --baseline) MIN_ALIAS="${2:-}"; shift 2 ;;
    --max) MAX_PARALLEL="${2:-}"; shift 2 ;;
    --runner-user) RUNNER_USER="${2:-}"; shift 2 ;;
    --base-dir) BASE_DIR="${2:-}"; shift 2 ;;
    --base-name) BASE_NAME="${2:-}"; shift 2 ;;
    --labels) LABELS="${2:-}"; shift 2 ;;
    --runner-version) RUNNER_VERSION="${2:-}"; shift 2 ;;
    --seed-dir) SEED_DIR="${2:-}"; shift 2 ;;
    --inventory-out) INVENTORY_OUT="${2:-}"; shift 2 ;;
    --status-out) STATUS_OUT="${2:-}"; shift 2 ;;
    --force-remove-busy) FORCE_REMOVE_BUSY="true"; shift 1 ;;
    --dry-run) DRY_RUN="true"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -n "$MIN_ALIAS" ]]; then
  if [[ "$MIN_POOL" != "1" && "$MIN_POOL" != "$MIN_ALIAS" ]]; then
    echo "--min and --baseline conflict; provide only one value." >&2
    exit 1
  fi
  MIN_POOL="$MIN_ALIAS"
fi

if [[ -z "$REPO" ]]; then
  echo "--repo is required." >&2
  usage
  exit 1
fi
if ! [[ "$MIN_POOL" =~ ^[0-9]+$ ]]; then
  echo "--min must be a non-negative integer." >&2
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
if [[ "$MIN_POOL" -gt "$MAX_PARALLEL" ]]; then
  echo "--min must be <= --max." >&2
  exit 1
fi
if [[ -n "$TARGET" && ! "$TARGET" =~ ^[0-9]+$ ]]; then
  echo "--target must be a non-negative integer when provided." >&2
  exit 1
fi

requested_target="${TARGET:-$MIN_POOL}"
effective_target="$requested_target"
if [[ "$effective_target" -lt "$MIN_POOL" ]]; then
  effective_target="$MIN_POOL"
fi
if [[ "$effective_target" -gt "$MAX_PARALLEL" ]]; then
  effective_target="$MAX_PARALLEL"
fi

require_cmd gh
require_cmd jq
require_cmd systemctl
require_cmd rsync

if [[ "${EUID}" -ne 0 ]]; then
  require_cmd sudo
  if [[ "$DRY_RUN" != "true" ]] && ! sudo -n true >/dev/null 2>&1; then
    echo "sudo auth is required. Run: sudo -v" >&2
    exit 1
  fi
else
  require_cmd runuser
fi

declare -a OPERATION_LINES=()
declare -a WARNING_LINES=()

echo "Scaling runner pool for ${REPO}"
echo "min=${MIN_POOL} max=${MAX_PARALLEL} requested_target=${requested_target} effective_target=${effective_target}"
echo "base_name=${BASE_NAME} base_dir=${BASE_DIR}"

# Build desired runner set.
declare -A DESIRED_MAP=()
declare -a DESIRED_NAMES=()
for ((i=1; i<=effective_target; i++)); do
  n="$(runner_name_for_index "$i")"
  DESIRED_MAP["$n"]=1
  DESIRED_NAMES+=("$n")
done

# Ensure desired runners exist and are running.
for n in "${DESIRED_NAMES[@]}"; do
  d="$(runner_dir_for_name "$n")"
  if [[ -d "$d" && -x "${d}/config.sh" && -f "${d}/.runner" ]]; then
    echo "exists (configured): ${n}"
    run_sudo bash -lc "cd '${d}' && ./svc.sh start || true"
    add_operation "ensure_started" "$n" "existing_config"
    continue
  fi
  if [[ -d "$d" && -x "${d}/config.sh" ]]; then
    echo "reconfigure (missing .runner): ${n}"
    configure_runner "$n" "reconfigure"
  else
    echo "create: ${n}"
    configure_runner "$n" "create"
  fi
done

# Scale down extra managed runners using deterministic order (highest index first).
mapfile -t managed_desc < <(managed_names_from_github_sorted desc)
for n in "${managed_desc[@]}"; do
  if [[ -z "${DESIRED_MAP[$n]+x}" ]]; then
    echo "remove: ${n}"
    remove_runner "$n"
  fi
done

# Inventory snapshots collected after apply phase.
mapfile -t GITHUB_MANAGED_NAMES < <(managed_names_from_github_sorted asc)
mapfile -t LOCAL_MANAGED_NAMES < <(managed_names_from_local_sorted)

declare -A GITHUB_MAP=()
declare -A LOCAL_MAP=()
for n in "${GITHUB_MANAGED_NAMES[@]}"; do GITHUB_MAP["$n"]=1; done
for n in "${LOCAL_MANAGED_NAMES[@]}"; do LOCAL_MAP["$n"]=1; done

if [[ -n "$INVENTORY_OUT" ]]; then
  emit_inventory
  echo "inventory_json=${INVENTORY_OUT}"
fi
if [[ -n "$STATUS_OUT" ]]; then
  emit_status
  echo "status_json=${STATUS_OUT}"
fi

echo "Done."
