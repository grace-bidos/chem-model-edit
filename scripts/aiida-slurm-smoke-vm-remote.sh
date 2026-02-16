#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

REMOTE_VM_HOST="${REMOTE_VM_HOST:-}"
REMOTE_VM_USER="${REMOTE_VM_USER:-}"
REMOTE_VM_PORT="${REMOTE_VM_PORT:-22}"
REMOTE_VM_REPO_DIR="${REMOTE_VM_REPO_DIR:-}"
REMOTE_VM_SMOKE_SCRIPT="${REMOTE_VM_SMOKE_SCRIPT:-scripts/aiida-slurm-smoke-vm.sh}"
REMOTE_VM_SSH_KEY="${REMOTE_VM_SSH_KEY:-}"
REMOTE_SSH_BIN="${REMOTE_SSH_BIN:-ssh}"
LOCAL_ARTIFACT_BASE="${LOCAL_ARTIFACT_BASE:-$ROOT_DIR/investigations/artifacts/gra-124}"
REMOTE_ARTIFACT_BASE="${REMOTE_ARTIFACT_BASE:-/tmp/gra-124-remote-vm-validation}"
RUN_ID="${GRA124_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
PROFILE="${AIIDA_PROFILE:-}"
COMPUTER_LABEL="${AIIDA_COMPUTER_LABEL:-}"
WORKDIR_TEMPLATE="${AIIDA_WORKDIR_TEMPLATE:-}"
DRY_RUN=0

declare -a SSH_OPTIONS=()
declare -a SSH_CMD=()

RUN_DIR=""
REMOTE_TARGET=""
REMOTE_RUN_NAME=""
REMOTE_ARTIFACT_DIR=""

usage() {
  cat <<'USAGE'
Usage: aiida-slurm-smoke-vm-remote.sh [options]

Run the VM Slurm+AiiDA smoke script remotely over SSH and collect evidence
artifacts to the local workspace.

Options:
  --host <hostname|user@hostname>   Remote VM host (env: REMOTE_VM_HOST)
  --user <username>                 SSH user (env: REMOTE_VM_USER)
  --port <port>                     SSH port (env: REMOTE_VM_PORT, default: 22)
  --repo-dir <path>                 Remote repository root (env: REMOTE_VM_REPO_DIR)
  --remote-script <path>            Remote smoke script path, absolute or repo-relative
                                    (env: REMOTE_VM_SMOKE_SCRIPT, default: scripts/aiida-slurm-smoke-vm.sh)
  --ssh-bin <path>                  SSH binary path (env: REMOTE_SSH_BIN, default: ssh)
  --ssh-key <path>                  SSH private key path (env: REMOTE_VM_SSH_KEY)
  --ssh-opt <value>                 Extra ssh -o option (repeatable)
  --profile <name>                  AiiDA profile passed to remote smoke script
  --computer-label <label>          AiiDA computer label passed to remote smoke script
  --workdir <path>                  AiiDA workdir template passed to remote smoke script
  --local-artifact-base <path>      Local artifact root
                                    (env: LOCAL_ARTIFACT_BASE, default: investigations/artifacts/gra-124)
  --remote-artifact-base <path>     Remote artifact root
                                    (env: REMOTE_ARTIFACT_BASE, default: /tmp/gra-124-remote-vm-validation)
  --run-id <id>                     Deterministic run id for local/remote artifact naming
                                    (env: GRA124_RUN_ID, default: UTC timestamp)
  --dry-run                         Print planned steps without executing SSH commands
  -h, --help                        Show this help

Exit codes:
  0  Remote smoke passed and artifacts were collected locally
  1  Local invocation/configuration error or hard remote script failure
  2  Blocked by remote prerequisite/connectivity/evidence collection issue
USAGE
}

log() {
  printf '%s\n' "$*"
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

have_command() {
  command -v "$1" >/dev/null 2>&1
}

to_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  elif [[ "$value" == "~/"* ]]; then
    printf '%s/%s\n' "$HOME" "${value#~/}"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

to_wsl_path_if_windows() {
  local value="$1"
  local drive=""
  local rest=""
  if [[ "$value" =~ ^([A-Za-z]):\\(.*)$ ]]; then
    drive="${BASH_REMATCH[1],,}"
    rest="${BASH_REMATCH[2]//\\//}"
    printf '/mnt/%s/%s\n' "$drive" "$rest"
    return 0
  fi
  if [[ "$value" =~ ^([A-Za-z]):/(.*)$ ]]; then
    drive="${BASH_REMATCH[1],,}"
    rest="${BASH_REMATCH[2]}"
    printf '/mnt/%s/%s\n' "$drive" "$rest"
    return 0
  fi
  printf '%s\n' "$value"
}

slugify() {
  local value="$1"
  local slug=""
  slug="$(printf '%s' "$value" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9._-]+/-/g; s/^-+//; s/-+$//')"
  if [[ -z "$slug" ]]; then
    slug="remote-vm"
  fi
  printf '%s\n' "$slug"
}

q() {
  printf '%q' "$1"
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --host)
        shift
        [[ $# -gt 0 ]] || die "--host requires a value"
        REMOTE_VM_HOST="$1"
        ;;
      --user)
        shift
        [[ $# -gt 0 ]] || die "--user requires a value"
        REMOTE_VM_USER="$1"
        ;;
      --port)
        shift
        [[ $# -gt 0 ]] || die "--port requires a value"
        REMOTE_VM_PORT="$1"
        ;;
      --repo-dir)
        shift
        [[ $# -gt 0 ]] || die "--repo-dir requires a value"
        REMOTE_VM_REPO_DIR="$1"
        ;;
      --remote-script)
        shift
        [[ $# -gt 0 ]] || die "--remote-script requires a value"
        REMOTE_VM_SMOKE_SCRIPT="$1"
        ;;
      --ssh-bin)
        shift
        [[ $# -gt 0 ]] || die "--ssh-bin requires a value"
        REMOTE_SSH_BIN="$1"
        ;;
      --ssh-key)
        shift
        [[ $# -gt 0 ]] || die "--ssh-key requires a value"
        REMOTE_VM_SSH_KEY="$1"
        ;;
      --ssh-opt)
        shift
        [[ $# -gt 0 ]] || die "--ssh-opt requires a value"
        SSH_OPTIONS+=("$1")
        ;;
      --profile)
        shift
        [[ $# -gt 0 ]] || die "--profile requires a value"
        PROFILE="$1"
        ;;
      --computer-label)
        shift
        [[ $# -gt 0 ]] || die "--computer-label requires a value"
        COMPUTER_LABEL="$1"
        ;;
      --workdir)
        shift
        [[ $# -gt 0 ]] || die "--workdir requires a value"
        WORKDIR_TEMPLATE="$1"
        ;;
      --local-artifact-base)
        shift
        [[ $# -gt 0 ]] || die "--local-artifact-base requires a value"
        LOCAL_ARTIFACT_BASE="$1"
        ;;
      --remote-artifact-base)
        shift
        [[ $# -gt 0 ]] || die "--remote-artifact-base requires a value"
        REMOTE_ARTIFACT_BASE="$1"
        ;;
      --run-id)
        shift
        [[ $# -gt 0 ]] || die "--run-id requires a value"
        RUN_ID="$1"
        ;;
      --dry-run)
        DRY_RUN=1
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown option: $1"
        ;;
    esac
    shift
  done

  [[ -n "$REMOTE_VM_HOST" ]] || die "remote host is required (--host or REMOTE_VM_HOST)"
  [[ -n "$REMOTE_VM_REPO_DIR" ]] || die "remote repo dir is required (--repo-dir or REMOTE_VM_REPO_DIR)"
  [[ "$REMOTE_VM_PORT" =~ ^[0-9]+$ ]] || die "remote port must be numeric: $REMOTE_VM_PORT"
  [[ -n "$RUN_ID" ]] || die "run id must not be empty"

  LOCAL_ARTIFACT_BASE="$(to_abs_path "$LOCAL_ARTIFACT_BASE")"
  REMOTE_VM_REPO_DIR="$(to_wsl_path_if_windows "$REMOTE_VM_REPO_DIR")"
  REMOTE_VM_REPO_DIR="${REMOTE_VM_REPO_DIR/#\~/$HOME}"
  REMOTE_SSH_BIN="$(to_wsl_path_if_windows "$REMOTE_SSH_BIN")"

  if [[ -n "$REMOTE_VM_SSH_KEY" ]]; then
    REMOTE_VM_SSH_KEY="$(to_wsl_path_if_windows "$REMOTE_VM_SSH_KEY")"
    if [[ "$REMOTE_VM_SSH_KEY" != /* && "$REMOTE_VM_SSH_KEY" != "~/"* ]]; then
      REMOTE_VM_SSH_KEY="$(to_abs_path "$REMOTE_VM_SSH_KEY")"
    fi
    REMOTE_VM_SSH_KEY="${REMOTE_VM_SSH_KEY/#\~/$HOME}"
  fi
}

build_remote_target() {
  if [[ "$REMOTE_VM_HOST" == *"@"* ]]; then
    REMOTE_TARGET="$REMOTE_VM_HOST"
    if [[ -n "$REMOTE_VM_USER" ]]; then
      log "note: --host already includes user; ignoring --user"
    fi
    return 0
  fi

  if [[ -n "$REMOTE_VM_USER" ]]; then
    REMOTE_TARGET="${REMOTE_VM_USER}@${REMOTE_VM_HOST}"
  else
    REMOTE_TARGET="$REMOTE_VM_HOST"
  fi
}

init_run_dir() {
  local target_slug=""
  target_slug="$(slugify "$REMOTE_TARGET")"
  REMOTE_RUN_NAME="${RUN_ID}-${target_slug}-slurm-aiida-smoke"
  REMOTE_ARTIFACT_DIR="${REMOTE_ARTIFACT_BASE%/}/${REMOTE_RUN_NAME}"
  RUN_DIR="${LOCAL_ARTIFACT_BASE%/}/${target_slug}/${REMOTE_RUN_NAME}"
  mkdir -p "$RUN_DIR"
}

write_metadata() {
  cat >"$RUN_DIR/run-metadata.txt" <<EOF
run_id=$RUN_ID
remote_target=$REMOTE_TARGET
remote_port=$REMOTE_VM_PORT
remote_repo_dir=$REMOTE_VM_REPO_DIR
remote_script=$REMOTE_VM_SMOKE_SCRIPT
remote_artifact_dir=$REMOTE_ARTIFACT_DIR
local_run_dir=$RUN_DIR
profile=${PROFILE:-<default>}
computer_label=${COMPUTER_LABEL:-<default>}
workdir_template=${WORKDIR_TEMPLATE:-<default>}
dry_run=$DRY_RUN
generated_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
}

build_ssh_cmd() {
  SSH_CMD=("$REMOTE_SSH_BIN" "-o" "BatchMode=yes" "-o" "ConnectTimeout=15" "-p" "$REMOTE_VM_PORT")
  if [[ -n "$REMOTE_VM_SSH_KEY" ]]; then
    SSH_CMD+=("-i" "$REMOTE_VM_SSH_KEY")
  fi
  local opt=""
  for opt in "${SSH_OPTIONS[@]}"; do
    SSH_CMD+=("-o" "$opt")
  done
}

run_capture() {
  local name="$1"
  shift
  local log_path="$RUN_DIR/${name}.log"

  log ""
  log "[run] $*"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run] %s\n' "$*" | tee "$log_path"
    return 0
  fi

  set +e
  "$@" 2>&1 | tee "$log_path"
  local status=${PIPESTATUS[0]}
  set -e

  if [[ "$status" -eq 0 ]]; then
    log "[ok] ${name}"
  else
    log "[fail] ${name} (exit=${status})"
  fi
  return "$status"
}

build_remote_smoke_cmd() {
  local remote_script_path="$REMOTE_VM_SMOKE_SCRIPT"
  if [[ "$remote_script_path" != /* ]]; then
    remote_script_path="${REMOTE_VM_REPO_DIR%/}/$remote_script_path"
  fi

  local -a smoke_cmd=(bash "$remote_script_path" --artifact-dir "$REMOTE_ARTIFACT_DIR")
  if [[ -n "$PROFILE" ]]; then
    smoke_cmd+=(--profile "$PROFILE")
  fi
  if [[ -n "$COMPUTER_LABEL" ]]; then
    smoke_cmd+=(--computer-label "$COMPUTER_LABEL")
  fi
  if [[ -n "$WORKDIR_TEMPLATE" ]]; then
    smoke_cmd+=(--workdir "$WORKDIR_TEMPLATE")
  fi

  printf 'set -euo pipefail; mkdir -p %s; cd %s; %s' \
    "$(q "$REMOTE_ARTIFACT_BASE")" \
    "$(q "$REMOTE_VM_REPO_DIR")" \
    "$(printf '%q ' "${smoke_cmd[@]}")"
}

collect_remote_artifacts() {
  local bundle="$RUN_DIR/remote-artifacts.tar.gz"
  local stderr_log="$RUN_DIR/collect-remote-artifacts.stderr.log"
  local extract_dir="$RUN_DIR/remote-artifacts"
  local tar_cmd=""

  tar_cmd="$(printf 'set -euo pipefail; if [ -d %s ]; then tar -C %s -czf - %s; else exit 3; fi' \
    "$(q "$REMOTE_ARTIFACT_DIR")" \
    "$(q "$REMOTE_ARTIFACT_BASE")" \
    "$(q "$REMOTE_RUN_NAME")")"

  log ""
  log "[run] collect remote artifacts -> $bundle"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run] %s %s "%s"\n' "${SSH_CMD[*]}" "$REMOTE_TARGET" "$tar_cmd" | tee "$RUN_DIR/collect-remote-artifacts.log"
    return 0
  fi

  set +e
  "${SSH_CMD[@]}" "$REMOTE_TARGET" "$tar_cmd" >"$bundle" 2>"$stderr_log"
  local status=$?
  set -e

  if [[ "$status" -ne 0 ]]; then
    log "[fail] collect remote artifacts (exit=${status})"
    return "$status"
  fi

  mkdir -p "$extract_dir"
  tar -xzf "$bundle" -C "$extract_dir"
  log "[ok] collect remote artifacts"
  return 0
}

write_blocker_guidance() {
  local reason="$1"
  cat >"$RUN_DIR/fallback-next-steps.txt" <<EOF
Blocked on remote Slurm+AiiDA VM smoke rerun.
reason=$reason

Suggested next checks:
1. SSH path/connectivity
   - ${REMOTE_SSH_BIN} -p ${REMOTE_VM_PORT} ${REMOTE_TARGET} "hostname"
2. Remote script path and repository root
   - test -x ${REMOTE_VM_REPO_DIR%/}/${REMOTE_VM_SMOKE_SCRIPT}
3. Run the smoke script directly on the VM for immediate diagnostics
   - cd ${REMOTE_VM_REPO_DIR}
   - ${REMOTE_VM_SMOKE_SCRIPT} --artifact-dir ${REMOTE_ARTIFACT_DIR}

Local evidence:
- run metadata: ${RUN_DIR}/run-metadata.txt
- remote preflight log: ${RUN_DIR}/remote-preflight.log
- remote smoke log: ${RUN_DIR}/remote-smoke.log
- artifact collection stderr (if any): ${RUN_DIR}/collect-remote-artifacts.stderr.log
EOF
}

validate_local_prerequisites() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    return 0
  fi
  [[ -x "$REMOTE_SSH_BIN" ]] || have_command "$REMOTE_SSH_BIN" || die "ssh binary not found: $REMOTE_SSH_BIN"
  have_command tar || die "tar is required for artifact collection"
}

main() {
  parse_args "$@"
  build_remote_target
  init_run_dir
  write_metadata
  validate_local_prerequisites
  build_ssh_cmd

  log "Remote VM validation rerun start"
  log "  run dir: $RUN_DIR"
  log "  remote target: $REMOTE_TARGET"
  log "  remote artifact dir: $REMOTE_ARTIFACT_DIR"

  local preflight_cmd=""
  preflight_cmd='set -euo pipefail; hostname; uname -a; whoami; command -v verdi || true; command -v scontrol || true; command -v sinfo || true'

  if ! run_capture "remote-preflight" "${SSH_CMD[@]}" "$REMOTE_TARGET" "$preflight_cmd"; then
    log "blocker: unable to reach remote VM over SSH preflight."
    write_blocker_guidance "ssh-preflight-failed"
    exit 2
  fi

  local remote_smoke_cmd=""
  remote_smoke_cmd="$(build_remote_smoke_cmd)"
  local smoke_status=0
  if run_capture "remote-smoke" "${SSH_CMD[@]}" "$REMOTE_TARGET" "$remote_smoke_cmd"; then
    smoke_status=0
  else
    smoke_status=$?
  fi

  local collect_status=0
  if collect_remote_artifacts; then
    collect_status=0
  else
    collect_status=$?
  fi

  local final_exit=0
  if [[ "$smoke_status" -eq 0 && "$collect_status" -eq 0 ]]; then
    final_exit=0
  elif [[ "$smoke_status" -eq 2 ]]; then
    final_exit=2
  elif [[ "$collect_status" -ne 0 ]]; then
    final_exit=2
  elif [[ "$smoke_status" -eq 1 ]]; then
    final_exit=1
  else
    final_exit=2
  fi

  if [[ "$final_exit" -ne 0 ]]; then
    if [[ "$collect_status" -ne 0 ]]; then
      log "blocker: failed to collect remote evidence artifacts."
      write_blocker_guidance "artifact-collection-failed"
    elif [[ "$smoke_status" -eq 2 ]]; then
      log "blocker: remote smoke reported runtime/prerequisite blockers."
      write_blocker_guidance "remote-smoke-blocked"
    else
      log "error: remote smoke exited with hard failure."
      write_blocker_guidance "remote-smoke-hard-failure"
    fi
  fi

  log ""
  log "Remote rerun summary"
  log "  run dir: $RUN_DIR"
  log "  remote target: $REMOTE_TARGET"
  log "  remote smoke exit: $smoke_status"
  log "  artifact collection exit: $collect_status"
  if [[ "$final_exit" -eq 0 ]]; then
    log "  result: success"
  elif [[ "$final_exit" -eq 2 ]]; then
    log "  result: blocked (see fallback-next-steps.txt)"
  else
    log "  result: hard failure"
  fi

  exit "$final_exit"
}

main "$@"
