#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

REMOTE_VM_HOST="${REMOTE_VM_HOST:-}"
REMOTE_VM_USER="${REMOTE_VM_USER:-}"
REMOTE_VM_PORT="${REMOTE_VM_PORT:-22}"
REMOTE_VM_REPO_DIR="${REMOTE_VM_REPO_DIR:-}"
REMOTE_VM_BOOTSTRAP_SCRIPT="${REMOTE_VM_BOOTSTRAP_SCRIPT:-scripts/aiida-vm-bootstrap.sh}"
REMOTE_VM_SMOKE_SCRIPT="${REMOTE_VM_SMOKE_SCRIPT:-scripts/aiida-slurm-smoke-vm.sh}"
REMOTE_AIIDA_ENV_FILE="${REMOTE_AIIDA_ENV_FILE:-ops/aiida-vm/aiida-vm.env}"
REMOTE_VM_SSH_KEY="${REMOTE_VM_SSH_KEY:-}"
REMOTE_SSH_BIN="${REMOTE_SSH_BIN:-ssh}"
LOCAL_ARTIFACT_BASE="${LOCAL_ARTIFACT_BASE:-$ROOT_DIR/investigations/artifacts/gra-126}"
REMOTE_ARTIFACT_BASE="${REMOTE_ARTIFACT_BASE:-/tmp/gra-126-promotion-run-vm}"
RUN_ID="${GRA126_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
PROFILE="${AIIDA_PROFILE:-}"
COMPUTER_LABEL="${AIIDA_COMPUTER_LABEL:-}"
WORKDIR_TEMPLATE="${AIIDA_WORKDIR_TEMPLATE:-}"
BOOTSTRAP_APPLY=0
BOOTSTRAP_INIT_PROFILE=0
DRY_RUN=0

declare -a SSH_OPTIONS=()
declare -a SSH_CMD=()

RUN_DIR=""
REMOTE_TARGET=""
REMOTE_RUN_NAME=""
REMOTE_GATE_DIR=""

usage() {
  cat <<'USAGE'
Usage: aiida-promotion-gate-vm-remote.sh [options]

Run promotion gate VM checks (gate 3/4 from docs/process/aiida-promotion-gate-vm.md)
on a remote VM over SSH and collect evidence artifacts locally.

Options:
  --host <hostname|user@hostname>   Remote VM host (env: REMOTE_VM_HOST)
  --user <username>                 SSH user (env: REMOTE_VM_USER)
  --port <port>                     SSH port (env: REMOTE_VM_PORT, default: 22)
  --repo-dir <path>                 Remote repository root (env: REMOTE_VM_REPO_DIR)
  --bootstrap-script <path>         Remote bootstrap script path, absolute or repo-relative
                                    (env: REMOTE_VM_BOOTSTRAP_SCRIPT, default: scripts/aiida-vm-bootstrap.sh)
  --smoke-script <path>             Remote slurm smoke script path, absolute or repo-relative
                                    (env: REMOTE_VM_SMOKE_SCRIPT, default: scripts/aiida-slurm-smoke-vm.sh)
  --env-file <path>                 Remote AiiDA VM env file path passed to bootstrap script
                                    (env: REMOTE_AIIDA_ENV_FILE, default: ops/aiida-vm/aiida-vm.env)
  --bootstrap-apply                 Add --apply to remote bootstrap command
  --bootstrap-init-profile          Add --init-profile to remote bootstrap command
  --ssh-bin <path>                  SSH binary path (env: REMOTE_SSH_BIN, default: ssh)
  --ssh-key <path>                  SSH private key path (env: REMOTE_VM_SSH_KEY)
  --ssh-opt <value>                 Extra ssh -o option (repeatable)
  --profile <name>                  AiiDA profile passed to remote smoke script
  --computer-label <label>          AiiDA computer label passed to remote smoke script
  --workdir <path>                  AiiDA workdir template passed to remote smoke script
  --local-artifact-base <path>      Local artifact root
                                    (env: LOCAL_ARTIFACT_BASE, default: investigations/artifacts/gra-126)
  --remote-artifact-base <path>     Remote artifact root
                                    (env: REMOTE_ARTIFACT_BASE, default: /tmp/gra-126-promotion-run-vm)
  --run-id <id>                     Deterministic run id for local/remote artifact naming
                                    (env: GRA126_RUN_ID, default: UTC timestamp)
  --dry-run                         Print planned steps without executing SSH commands
  -h, --help                        Show this help

Exit codes:
  0  Promotion VM checks passed and artifacts were collected locally
  1  Local invocation/configuration error or hard remote script failure
  2  Blocked by remote prerequisite/connectivity/runtime/evidence collection issue
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

to_windows_path_if_wsl() {
  local value="$1"
  local drive=""
  local rest=""
  if [[ "$value" =~ ^/mnt/([A-Za-z])/(.*)$ ]]; then
    drive="${BASH_REMATCH[1]^}"
    rest="${BASH_REMATCH[2]}"
    printf '%s:/%s\n' "$drive" "$rest"
    return 0
  fi
  printf '%s\n' "$value"
}

is_windows_ssh_bin() {
  local value
  value="$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')"
  [[ "$value" == *"ssh.exe" ]]
}

normalize_ssh_key_path_for_ssh_bin() {
  local key_path="$1"
  local ssh_bin_path="$2"
  if is_windows_ssh_bin "$ssh_bin_path"; then
    printf '%s\n' "$(to_windows_path_if_wsl "$key_path")"
    return 0
  fi
  printf '%s\n' "$key_path"
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
      --bootstrap-script)
        shift
        [[ $# -gt 0 ]] || die "--bootstrap-script requires a value"
        REMOTE_VM_BOOTSTRAP_SCRIPT="$1"
        ;;
      --smoke-script)
        shift
        [[ $# -gt 0 ]] || die "--smoke-script requires a value"
        REMOTE_VM_SMOKE_SCRIPT="$1"
        ;;
      --env-file)
        shift
        [[ $# -gt 0 ]] || die "--env-file requires a value"
        REMOTE_AIIDA_ENV_FILE="$1"
        ;;
      --bootstrap-apply)
        BOOTSTRAP_APPLY=1
        ;;
      --bootstrap-init-profile)
        BOOTSTRAP_INIT_PROFILE=1
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
    REMOTE_VM_SSH_KEY="$(normalize_ssh_key_path_for_ssh_bin "$REMOTE_VM_SSH_KEY" "$REMOTE_SSH_BIN")"
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
  REMOTE_RUN_NAME="${RUN_ID}-${target_slug}-promotion-gate-vm"
  REMOTE_GATE_DIR="${REMOTE_ARTIFACT_BASE%/}/${REMOTE_RUN_NAME}"
  RUN_DIR="${LOCAL_ARTIFACT_BASE%/}/${target_slug}/${REMOTE_RUN_NAME}"
  mkdir -p "$RUN_DIR"
}

write_metadata() {
  cat >"$RUN_DIR/run-metadata.txt" <<EOF_META
run_id=$RUN_ID
remote_target=$REMOTE_TARGET
remote_port=$REMOTE_VM_PORT
remote_repo_dir=$REMOTE_VM_REPO_DIR
bootstrap_script=$REMOTE_VM_BOOTSTRAP_SCRIPT
smoke_script=$REMOTE_VM_SMOKE_SCRIPT
remote_env_file=$REMOTE_AIIDA_ENV_FILE
remote_gate_dir=$REMOTE_GATE_DIR
local_run_dir=$RUN_DIR
profile=${PROFILE:-<default>}
computer_label=${COMPUTER_LABEL:-<default>}
workdir_template=${WORKDIR_TEMPLATE:-<default>}
bootstrap_apply=$BOOTSTRAP_APPLY
bootstrap_init_profile=$BOOTSTRAP_INIT_PROFILE
dry_run=$DRY_RUN
generated_at_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF_META
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

build_remote_promotion_cmd() {
  local bootstrap_script_path="$REMOTE_VM_BOOTSTRAP_SCRIPT"
  local smoke_script_path="$REMOTE_VM_SMOKE_SCRIPT"

  if [[ "$bootstrap_script_path" != /* ]]; then
    bootstrap_script_path="${REMOTE_VM_REPO_DIR%/}/$bootstrap_script_path"
  fi
  if [[ "$smoke_script_path" != /* ]]; then
    smoke_script_path="${REMOTE_VM_REPO_DIR%/}/$smoke_script_path"
  fi

  local -a bootstrap_cmd=(bash "$bootstrap_script_path" --env-file "$REMOTE_AIIDA_ENV_FILE" --sanity-check)
  if [[ "$BOOTSTRAP_APPLY" -eq 1 ]]; then
    bootstrap_cmd=(bash "$bootstrap_script_path" --apply --env-file "$REMOTE_AIIDA_ENV_FILE" --sanity-check)
  fi
  if [[ "$BOOTSTRAP_INIT_PROFILE" -eq 1 ]]; then
    bootstrap_cmd+=(--init-profile)
  fi

  local -a smoke_cmd=(bash "$smoke_script_path" --artifact-dir "$REMOTE_GATE_DIR/slurm-smoke")
  if [[ -n "$PROFILE" ]]; then
    smoke_cmd+=(--profile "$PROFILE")
  fi
  if [[ -n "$COMPUTER_LABEL" ]]; then
    smoke_cmd+=(--computer-label "$COMPUTER_LABEL")
  fi
  if [[ -n "$WORKDIR_TEMPLATE" ]]; then
    smoke_cmd+=(--workdir "$WORKDIR_TEMPLATE")
  fi

  local bootstrap_q=""
  local smoke_q=""
  bootstrap_q="$(printf '%q ' "${bootstrap_cmd[@]}")"
  smoke_q="$(printf '%q ' "${smoke_cmd[@]}")"

  cat <<REMOTE_CMD
set -euo pipefail
mkdir -p $(q "$REMOTE_GATE_DIR")
cd $(q "$REMOTE_VM_REPO_DIR")

set +e
${bootstrap_q}2>&1 | tee $(q "$REMOTE_GATE_DIR/03-vm-bootstrap.log")
bootstrap_status=\${PIPESTATUS[0]}

if command -v verdi >/dev/null 2>&1; then
  verdi profile list 2>&1 | tee $(q "$REMOTE_GATE_DIR/03-verdi-profile-list.log")
  profile_status=\${PIPESTATUS[0]}
  verdi status 2>&1 | tee $(q "$REMOTE_GATE_DIR/03-verdi-status.log")
  verdi_status=\${PIPESTATUS[0]}
else
  printf 'verdi command missing on remote VM\n' | tee $(q "$REMOTE_GATE_DIR/03-verdi-profile-list.log")
  printf 'verdi command missing on remote VM\n' | tee $(q "$REMOTE_GATE_DIR/03-verdi-status.log")
  profile_status=127
  verdi_status=127
fi

if [ "\$bootstrap_status" -eq 0 ]; then
  ${smoke_q}2>&1 | tee $(q "$REMOTE_GATE_DIR/04-slurm-smoke-driver.log")
  smoke_status=\${PIPESTATUS[0]}
else
  printf 'skip slurm smoke because bootstrap failed (exit=%s)\n' "\$bootstrap_status" \
    | tee $(q "$REMOTE_GATE_DIR/04-slurm-smoke-driver.log")
  smoke_status=3
fi

cat > $(q "$REMOTE_GATE_DIR/gate-status.env") <<STATUS_EOF
bootstrap_exit=\$bootstrap_status
verdi_profile_list_exit=\$profile_status
verdi_status_exit=\$verdi_status
slurm_smoke_exit=\$smoke_status
STATUS_EOF

cat > $(q "$REMOTE_GATE_DIR/04-slurm-smoke-summary.md") <<SUMMARY_EOF
# Gate 4 Slurm Smoke Summary

- driver_exit: \$smoke_status
- evidence_dir: $REMOTE_GATE_DIR/slurm-smoke
- interpretation:
  - \`0\`: pass
  - \`2\`: blocked (runtime/prerequisite issue)
  - \`1\`: hard failure
  - \`3\`: skipped because gate 3 failed
SUMMARY_EOF

set -e
if [ "\$bootstrap_status" -ne 0 ]; then
  exit 20
fi
if [ "\$profile_status" -ne 0 ] || [ "\$verdi_status" -ne 0 ]; then
  exit 21
fi
if [ "\$smoke_status" -eq 2 ]; then
  exit 22
fi
if [ "\$smoke_status" -ne 0 ]; then
  exit 23
fi
exit 0
REMOTE_CMD
}

collect_remote_artifacts() {
  local bundle="$RUN_DIR/remote-artifacts.tar.gz"
  local stderr_log="$RUN_DIR/collect-remote-artifacts.stderr.log"
  local extract_dir="$RUN_DIR/remote-artifacts"
  local tar_cmd=""

  tar_cmd="$(printf 'set -euo pipefail; if [ -d %s ]; then tar -C %s -czf - %s; else exit 3; fi' \
    "$(q "$REMOTE_GATE_DIR")" \
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

interpret_remote_status() {
  local status="$1"
  case "$status" in
    0) printf 'passed' ;;
    20) printf 'gate-3-bootstrap-failed' ;;
    21) printf 'gate-3-verdi-check-failed' ;;
    22) printf 'gate-4-slurm-blocked' ;;
    23) printf 'gate-4-slurm-hard-failure' ;;
    255) printf 'ssh-command-failed' ;;
    *) printf 'unexpected-remote-exit-%s' "$status" ;;
  esac
}

write_local_summary() {
  local remote_status="$1"
  local collect_status="$2"
  local final_exit="$3"
  local reason=""

  reason="$(interpret_remote_status "$remote_status")"

  cat >"$RUN_DIR/promotion-gate-vm-summary.md" <<EOF_SUM
# Promotion Gate VM Remote Summary

- run_id: $RUN_ID
- remote_target: $REMOTE_TARGET
- remote_status: $remote_status ($reason)
- artifact_collection_status: $collect_status
- final_exit: $final_exit
- local_run_dir: $RUN_DIR
- remote_gate_dir: $REMOTE_GATE_DIR

Gate interpretation:
- gate 3 passes only when \`03-vm-bootstrap.log\`, \`03-verdi-profile-list.log\`, and \`03-verdi-status.log\` are successful.
- gate 4 passes only when \`04-slurm-smoke-driver.log\` exits \`0\` and slurm smoke artifacts exist.
- \`final_exit=2\` means \`HOLD\` (blocked) for promotion decision.
- \`final_exit=1\` means hard-failure path requiring script/operator fix before rerun.
EOF_SUM
}

write_blocker_guidance() {
  local reason="$1"
  cat >"$RUN_DIR/fallback-next-steps.txt" <<EOF_BLOCK
Blocked on PostgreSQL+RabbitMQ promotion gate VM rerun.
reason=$reason

Suggested checks:
1. SSH path/connectivity
   - ${REMOTE_SSH_BIN} -p ${REMOTE_VM_PORT} ${REMOTE_TARGET} "hostname"
2. Remote script/env paths
   - test -f ${REMOTE_VM_REPO_DIR%/}/${REMOTE_VM_BOOTSTRAP_SCRIPT}
   - test -f ${REMOTE_VM_REPO_DIR%/}/${REMOTE_VM_SMOKE_SCRIPT}
   - test -f ${REMOTE_VM_REPO_DIR%/}/${REMOTE_AIIDA_ENV_FILE}
3. Run gate 3 directly on VM
   - cd ${REMOTE_VM_REPO_DIR}
   - ${REMOTE_VM_BOOTSTRAP_SCRIPT} --env-file ${REMOTE_AIIDA_ENV_FILE} --sanity-check
4. If gate 3 passes, rerun gate 4 directly on VM
   - ${REMOTE_VM_SMOKE_SCRIPT} --artifact-dir ${REMOTE_GATE_DIR}/slurm-smoke

Local evidence:
- run metadata: ${RUN_DIR}/run-metadata.txt
- remote preflight log: ${RUN_DIR}/remote-preflight.log
- remote promotion run log: ${RUN_DIR}/remote-promotion-run.log
- summary: ${RUN_DIR}/promotion-gate-vm-summary.md
- collected remote bundle: ${RUN_DIR}/remote-artifacts.tar.gz
EOF_BLOCK
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

  log "Promotion gate VM rerun start"
  log "  run dir: $RUN_DIR"
  log "  remote target: $REMOTE_TARGET"
  log "  remote gate dir: $REMOTE_GATE_DIR"

  local preflight_cmd=""
  preflight_cmd="$(printf 'set -euo pipefail; cd %s; test -f %s; test -f %s; test -f %s; hostname; uname -a; whoami; command -v verdi || true; command -v scontrol || true; command -v sinfo || true' \
    "$(q "$REMOTE_VM_REPO_DIR")" \
    "$(q "$REMOTE_VM_BOOTSTRAP_SCRIPT")" \
    "$(q "$REMOTE_VM_SMOKE_SCRIPT")" \
    "$(q "$REMOTE_AIIDA_ENV_FILE")")"

  if ! run_capture "remote-preflight" "${SSH_CMD[@]}" "$REMOTE_TARGET" "$preflight_cmd"; then
    log "blocker: unable to reach or preflight-check remote VM over SSH."
    write_blocker_guidance "ssh-or-preflight-failed"
    write_local_summary 255 0 2
    exit 2
  fi

  local remote_cmd=""
  remote_cmd="$(build_remote_promotion_cmd)"

  local remote_status=0
  if run_capture "remote-promotion-run" "${SSH_CMD[@]}" "$REMOTE_TARGET" "$remote_cmd"; then
    remote_status=0
  else
    remote_status=$?
  fi

  local collect_status=0
  if collect_remote_artifacts; then
    collect_status=0
  else
    collect_status=$?
  fi

  local final_exit=0
  if [[ "$remote_status" -eq 0 && "$collect_status" -eq 0 ]]; then
    final_exit=0
  elif [[ "$collect_status" -ne 0 ]]; then
    final_exit=2
  elif [[ "$remote_status" -eq 20 || "$remote_status" -eq 21 || "$remote_status" -eq 22 || "$remote_status" -eq 255 ]]; then
    final_exit=2
  elif [[ "$remote_status" -eq 23 ]]; then
    final_exit=1
  else
    final_exit=1
  fi

  write_local_summary "$remote_status" "$collect_status" "$final_exit"

  if [[ "$final_exit" -ne 0 ]]; then
    if [[ "$collect_status" -ne 0 ]]; then
      log "blocker: failed to collect remote evidence artifacts."
      write_blocker_guidance "artifact-collection-failed"
    elif [[ "$remote_status" -eq 20 || "$remote_status" -eq 21 ]]; then
      log "blocker: gate 3 (bootstrap/verdi checks) failed on remote VM."
      write_blocker_guidance "gate-3-failed"
    elif [[ "$remote_status" -eq 22 ]]; then
      log "blocker: gate 4 slurm smoke reported blocked state."
      write_blocker_guidance "gate-4-blocked"
    else
      log "error: remote promotion run exited with hard failure."
      write_blocker_guidance "remote-hard-failure"
    fi
  fi

  log ""
  log "Promotion gate VM rerun summary"
  log "  run dir: $RUN_DIR"
  log "  remote target: $REMOTE_TARGET"
  log "  remote status: $remote_status"
  log "  artifact collection status: $collect_status"
  if [[ "$final_exit" -eq 0 ]]; then
    log "  result: success"
  elif [[ "$final_exit" -eq 2 ]]; then
    log "  result: blocked (HOLD; see fallback-next-steps.txt)"
  else
    log "  result: hard failure"
  fi

  exit "$final_exit"
}

main "$@"
