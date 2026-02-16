#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MODE="vm"
SLURM_CONF=""
CGROUP_CONF=""
MUNGE_KEY=""
RUN_RUNTIME_CHECKS=1
CHECK_MUNGE_KEY_PERMS=1

declare -a FAILURES=()
declare -a BLOCKERS=()
declare -a WARNINGS=()

usage() {
  cat <<'USAGE'
Usage: validate-slurm-vm-control-plane.sh [options]

Validate Slurm VM control-plane prerequisites and config health.

Modes:
  vm       Validate a real VM host (default).
  offline  Deterministic static validation using repository fixtures.

Options:
  --mode <vm|offline>      Validation mode.
  --slurm-conf <path>      Path to slurm.conf.
  --cgroup-conf <path>     Path to cgroup.conf.
  --munge-key <path>       Path to munge.key.
  --skip-runtime           Do not run runtime checks.
  --skip-key-perms         Do not enforce munge key permissions.
  -h, --help               Show this help.

Exit codes:
  0  All checks passed.
  1  Deterministic validation failures were found.
  2  Runtime/prerequisite blockers were found (static checks still executed).
USAGE
}

log() {
  printf '%s\n' "$*"
}

add_failure() {
  FAILURES+=("$1")
}

add_blocker() {
  BLOCKERS+=("$1")
}

add_warning() {
  WARNINGS+=("$1")
}

have_command() {
  command -v "$1" >/dev/null 2>&1
}

to_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

set_defaults() {
  if [[ "$MODE" == "offline" ]]; then
    SLURM_CONF="${SLURM_CONF:-$ROOT_DIR/ops/slurm-vm/examples/slurm.conf}"
    CGROUP_CONF="${CGROUP_CONF:-$ROOT_DIR/ops/slurm-vm/examples/cgroup.conf}"
    MUNGE_KEY="${MUNGE_KEY:-$ROOT_DIR/ops/slurm-vm/examples/munge.key.example}"
    RUN_RUNTIME_CHECKS=0
    CHECK_MUNGE_KEY_PERMS=0
  else
    SLURM_CONF="${SLURM_CONF:-/etc/slurm/slurm.conf}"
    CGROUP_CONF="${CGROUP_CONF:-/etc/slurm/cgroup.conf}"
    MUNGE_KEY="${MUNGE_KEY:-/etc/munge/munge.key}"
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --mode)
        shift
        [[ $# -gt 0 ]] || { log "error: --mode requires a value"; exit 1; }
        MODE="$1"
        ;;
      --slurm-conf)
        shift
        [[ $# -gt 0 ]] || { log "error: --slurm-conf requires a value"; exit 1; }
        SLURM_CONF="$1"
        ;;
      --cgroup-conf)
        shift
        [[ $# -gt 0 ]] || { log "error: --cgroup-conf requires a value"; exit 1; }
        CGROUP_CONF="$1"
        ;;
      --munge-key)
        shift
        [[ $# -gt 0 ]] || { log "error: --munge-key requires a value"; exit 1; }
        MUNGE_KEY="$1"
        ;;
      --skip-runtime)
        RUN_RUNTIME_CHECKS=0
        ;;
      --skip-key-perms)
        CHECK_MUNGE_KEY_PERMS=0
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        log "error: unknown option: $1"
        usage
        exit 1
        ;;
    esac
    shift
  done

  if [[ "$MODE" != "vm" && "$MODE" != "offline" ]]; then
    log "error: unsupported --mode value: $MODE"
    exit 1
  fi

  set_defaults
  SLURM_CONF="$(to_abs_path "$SLURM_CONF")"
  CGROUP_CONF="$(to_abs_path "$CGROUP_CONF")"
  MUNGE_KEY="$(to_abs_path "$MUNGE_KEY")"
}

trim() {
  local value="$1"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  printf '%s\n' "$value"
}

conf_value() {
  local key="$1"
  local file="$2"
  awk -F '=' -v k="$key" '
    BEGIN { IGNORECASE = 0 }
    /^[[:space:]]*#/ { next }
    $1 ~ ("^[[:space:]]*" k "[[:space:]]*$") {
      val=$2
      sub(/^[[:space:]]+/, "", val)
      sub(/[[:space:]]+$/, "", val)
      print val
      exit
    }
  ' "$file"
}

check_commands() {
  local required=(munge unmunge slurmctld slurmd scontrol sinfo)
  for cmd in "${required[@]}"; do
    if ! have_command "$cmd"; then
      add_blocker "missing required command: $cmd"
    fi
  done

  if ! have_command systemctl; then
    add_warning "systemctl not found; service enable/active checks will be skipped"
  fi
}

check_paths() {
  [[ -f "$SLURM_CONF" ]] || add_blocker "missing slurm.conf: $SLURM_CONF"
  [[ -f "$CGROUP_CONF" ]] || add_blocker "missing cgroup.conf: $CGROUP_CONF"
  [[ -f "$MUNGE_KEY" ]] || add_blocker "missing munge key: $MUNGE_KEY"
}

validate_slurm_conf_static() {
  [[ -f "$SLURM_CONF" ]] || return 0

  local keys=(
    ClusterName
    SlurmctldHost
    SlurmctldPort
    SlurmUser
    StateSaveLocation
    SlurmdSpoolDir
    ProctrackType
    TaskPlugin
    AuthType
  )

  local key value
  for key in "${keys[@]}"; do
    value="$(trim "$(conf_value "$key" "$SLURM_CONF")")"
    if [[ -z "$value" ]]; then
      add_failure "slurm.conf missing required key: $key"
    fi
  done

  local auth_type
  auth_type="$(trim "$(conf_value "AuthType" "$SLURM_CONF")")"
  if [[ -n "$auth_type" && "$auth_type" != "auth/munge" ]]; then
    add_failure "slurm.conf AuthType must be auth/munge (found: $auth_type)"
  fi

  if ! grep -Eq '^[[:space:]]*NodeName[[:space:]]*=' "$SLURM_CONF"; then
    add_failure "slurm.conf must include at least one NodeName entry"
  fi

  if ! grep -Eq '^[[:space:]]*PartitionName[[:space:]]*=' "$SLURM_CONF"; then
    add_failure "slurm.conf must include at least one PartitionName entry"
  fi
}

validate_cgroup_conf_static() {
  [[ -f "$CGROUP_CONF" ]] || return 0

  local keys=(CgroupPlugin ConstrainCores ConstrainRAMSpace)
  local key value
  for key in "${keys[@]}"; do
    value="$(trim "$(conf_value "$key" "$CGROUP_CONF")")"
    if [[ -z "$value" ]]; then
      add_failure "cgroup.conf missing required key: $key"
    fi
  done
}

check_munge_key_permissions() {
  [[ "$CHECK_MUNGE_KEY_PERMS" -eq 1 ]] || return 0
  [[ -f "$MUNGE_KEY" ]] || return 0

  if ! have_command stat; then
    add_warning "stat not found; skipping munge key ownership/permission checks"
    return 0
  fi

  local mode owner group
  mode="$(stat -c '%a' "$MUNGE_KEY" 2>/dev/null || true)"
  owner="$(stat -c '%U' "$MUNGE_KEY" 2>/dev/null || true)"
  group="$(stat -c '%G' "$MUNGE_KEY" 2>/dev/null || true)"

  if [[ "$mode" != "400" ]]; then
    add_failure "munge key mode must be 400 (found: ${mode:-unknown})"
  fi
  if [[ "$owner" != "munge" || "$group" != "munge" ]]; then
    add_failure "munge key owner/group must be munge:munge (found: ${owner:-unknown}:${group:-unknown})"
  fi
}

run_runtime_check() {
  local label="$1"
  shift
  if "$@" >/tmp/validate-slurm-vm-control-plane.runtime.log 2>&1; then
    log "[ok] $label"
  else
    local output
    output="$(tr '\n' ' ' </tmp/validate-slurm-vm-control-plane.runtime.log | sed 's/[[:space:]]\+/ /g' | cut -c1-220)"
    add_blocker "$label failed: $output"
    log "[blocked] $label"
  fi
}

validate_runtime() {
  [[ "$RUN_RUNTIME_CHECKS" -eq 1 ]] || {
    log "runtime checks: skipped"
    return 0
  }

  local missing_runtime_tools=0
  local cmd
  for cmd in munge unmunge slurmctld slurmd scontrol; do
    if ! have_command "$cmd"; then
      missing_runtime_tools=1
    fi
  done
  if [[ "$missing_runtime_tools" -eq 1 ]]; then
    log "runtime checks: skipped due to missing commands"
    return 0
  fi

  if [[ ! -f "$SLURM_CONF" ]]; then
    log "runtime checks: skipped because slurm.conf is missing"
    return 0
  fi

  run_runtime_check "munge encode/decode" bash -lc "munge -n | unmunge >/dev/null"
  run_runtime_check "slurmctld config test" slurmctld -t -f "$SLURM_CONF"
  run_runtime_check "slurmd config test" slurmd -t -f "$SLURM_CONF"
  run_runtime_check "scontrol ping" scontrol ping

  if have_command systemctl; then
    local service
    for service in munge slurmctld slurmd; do
      if ! systemctl is-enabled --quiet "$service"; then
        add_blocker "service not enabled: $service"
      fi
      if ! systemctl is-active --quiet "$service"; then
        add_blocker "service not active: $service"
      fi
    done
  fi
}

report_results() {
  log ""
  log "mode: $MODE"
  log "slurm.conf: $SLURM_CONF"
  log "cgroup.conf: $CGROUP_CONF"
  log "munge key: $MUNGE_KEY"
  log ""

  if [[ "${#FAILURES[@]}" -gt 0 ]]; then
    log "deterministic validation failures:"
    printf '  - %s\n' "${FAILURES[@]}"
  else
    log "deterministic validation failures: none"
  fi

  if [[ "${#BLOCKERS[@]}" -gt 0 ]]; then
    log "runtime/prerequisite blockers:"
    printf '  - %s\n' "${BLOCKERS[@]}"
  else
    log "runtime/prerequisite blockers: none"
  fi

  if [[ "${#WARNINGS[@]}" -gt 0 ]]; then
    log "warnings:"
    printf '  - %s\n' "${WARNINGS[@]}"
  fi

  if [[ "${#FAILURES[@]}" -gt 0 ]]; then
    return 1
  fi
  if [[ "${#BLOCKERS[@]}" -gt 0 ]]; then
    return 2
  fi
  return 0
}

main() {
  parse_args "$@"
  check_commands
  check_paths
  validate_slurm_conf_static
  validate_cgroup_conf_static
  check_munge_key_permissions
  validate_runtime
  report_results
}

main "$@"
