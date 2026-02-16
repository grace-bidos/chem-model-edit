#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/investigations/artifacts/gra-89"
PROFILE="${AIIDA_PROFILE:-}"
COMPUTER_LABEL="gra89-vm-slurm"
WORKDIR_TEMPLATE="/tmp/aiida-gra89-slurm-{username}"
EXIT_CODE=0
RUN_DIR=""

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Run a VM-targeted AiiDA -> Slurm smoke execution.
This script assumes Slurm and AiiDA are already installed on the VM.

Options:
  --profile <name>        AiiDA profile name (default: current default profile)
  --computer-label <id>   AiiDA computer label (default: gra89-vm-slurm)
  --workdir <path>        Remote workdir template (default: /tmp/aiida-gra89-slurm-{username})
  --artifact-dir <path>   Output directory for logs (default: investigations/artifacts/gra-89)
  -h, --help              Show help

Exit codes:
  0  smoke passed
  1  hard failure (invalid arguments or unexpected script error)
  2  blocked (missing runtime/tooling/profile or Slurm control plane issue)
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
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
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
      --artifact-dir)
        shift
        [[ $# -gt 0 ]] || die "--artifact-dir requires a value"
        ARTIFACT_DIR="$1"
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

  ARTIFACT_DIR="$(to_abs_path "$ARTIFACT_DIR")"
}

run_capture() {
  local name="$1"
  shift
  local log_path="$RUN_DIR/${name}.log"

  log ""
  log "[run] $*"
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

computer_test_clean() {
  local log_path="$RUN_DIR/aiida-computer-test-slurm.log"
  if [[ ! -f "$log_path" ]]; then
    return 1
  fi

  if grep -Eq '\[Failed\]|tests failed' "$log_path"; then
    return 1
  fi
  return 0
}

resolve_profile() {
  if [[ -n "$PROFILE" ]]; then
    return 0
  fi

  set +e
  PROFILE="$(verdi profile show 2>/dev/null | tr -d '[:space:]')"
  local status=$?
  set -e

  if [[ "$status" -ne 0 || -z "$PROFILE" ]]; then
    return 1
  fi
  return 0
}

profile_exists() {
  verdi profile show "$PROFILE" >/dev/null 2>&1
}

computer_exists() {
  verdi -p "$PROFILE" computer show "$COMPUTER_LABEL" >/dev/null 2>&1
}

ensure_slurm_computer() {
  if computer_exists; then
    log "slurm computer already exists: ${COMPUTER_LABEL}"
    return 0
  fi

  run_capture "aiida-computer-setup-slurm" \
    verdi -p "$PROFILE" computer setup --non-interactive \
    --label "$COMPUTER_LABEL" \
    --hostname localhost \
    --transport core.local \
    --scheduler core.slurm \
    --work-dir "$WORKDIR_TEMPLATE"

  run_capture "aiida-computer-configure-local" \
    verdi -p "$PROFILE" computer configure core.local --non-interactive \
    --safe-interval 0 "$COMPUTER_LABEL"
}

collect_blocker_evidence() {
  log ""
  log "Collecting fallback diagnostics..."

  run_capture "diag-which-verdi" bash -lc 'command -v verdi'
  run_capture "diag-verdi-version" verdi --version
  run_capture "diag-verdi-profile-list" verdi profile list
  run_capture "diag-slurm-scontrol-ping" scontrol ping
  run_capture "diag-slurm-sinfo" sinfo
  run_capture "diag-slurm-squeue-user" squeue -u "$(id -un)"
  run_capture "diag-slurm-sbatch-version" sbatch --version

  cat >"$RUN_DIR/fallback-next-steps.txt" <<TXT
Blocked on AiiDA -> Slurm smoke.

Suggested next checks on the VM:
1. Verify AiiDA installation and active profile:
   - command -v verdi
   - verdi --version
   - verdi profile list
2. Verify Slurm controller reachability:
   - scontrol ping
   - sinfo
3. If Slurm is reachable but AiiDA test still fails:
   - verdi -p ${PROFILE:-<profile>} computer test ${COMPUTER_LABEL}
   - inspect this run directory logs for scheduler/transport details.
TXT

  log "fallback guidance: $RUN_DIR/fallback-next-steps.txt"
}

main() {
  parse_args "$@"

  mkdir -p "$ARTIFACT_DIR"
  RUN_DIR="$ARTIFACT_DIR/$(date +%Y%m%d-%H%M%S)-vm-slurm"
  mkdir -p "$RUN_DIR"

  log "AiiDA -> Slurm VM smoke start"
  log "  run dir: ${RUN_DIR}"

  run_capture "env-hostname" hostname || true
  run_capture "env-uname" uname -a || true
  run_capture "env-whoami" whoami || true

  local missing=()
  for cmd in verdi sinfo sbatch squeue scontrol; do
    if ! have_command "$cmd"; then
      missing+=("$cmd")
    fi
  done

  if [[ ${#missing[@]} -gt 0 ]]; then
    log "blocker: missing required commands: ${missing[*]}"
    EXIT_CODE=2
    collect_blocker_evidence || true
    log ""
    log "Smoke summary"
    log "  artifacts: ${RUN_DIR}"
    log "  result: blocked (missing runtime prerequisites)"
    exit "$EXIT_CODE"
  fi

  run_capture "env-verdi-version" verdi --version || true
  run_capture "env-slurm-version" sbatch --version || true

  if ! resolve_profile; then
    log "blocker: unable to resolve an AiiDA profile."
    EXIT_CODE=2
    collect_blocker_evidence || true
    exit "$EXIT_CODE"
  fi

  log "  profile: ${PROFILE}"
  log "  computer label: ${COMPUTER_LABEL}"

  if ! profile_exists; then
    log "blocker: AiiDA profile not found: ${PROFILE}"
    EXIT_CODE=2
    collect_blocker_evidence || true
    exit "$EXIT_CODE"
  fi

  if ! run_capture "slurm-scontrol-ping" scontrol ping; then
    log "blocker: slurm control plane check failed (scontrol ping)."
    EXIT_CODE=2
  fi

  if ! run_capture "slurm-sinfo" sinfo; then
    log "blocker: slurm scheduler listing failed (sinfo)."
    EXIT_CODE=2
  fi

  ensure_slurm_computer

  if run_capture "aiida-computer-test-slurm" verdi -p "$PROFILE" computer test "$COMPUTER_LABEL" && computer_test_clean; then
    EXIT_CODE=0
  else
    log "blocker: verdi computer test failed for ${COMPUTER_LABEL}."
    EXIT_CODE=2
  fi

  if [[ "$EXIT_CODE" -eq 2 ]]; then
    collect_blocker_evidence || true
  fi

  log ""
  log "Smoke summary"
  log "  profile: ${PROFILE}"
  log "  computer: ${COMPUTER_LABEL}"
  log "  artifacts: ${RUN_DIR}"

  if [[ "$EXIT_CODE" -eq 0 ]]; then
    log "  result: success"
  else
    log "  result: blocked (see logs + fallback-next-steps.txt)"
  fi

  exit "$EXIT_CODE"
}

main "$@"
