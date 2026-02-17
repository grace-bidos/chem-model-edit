#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ARTIFACT_BASE="${TIER2_ARTIFACT_BASE:-$ROOT_DIR/investigations/artifacts/gra-130}"
RUN_ID="${GRA130_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
VALIDATE_SCRIPT="${TIER2_VALIDATE_SCRIPT:-$ROOT_DIR/scripts/validate-slurm-vm-control-plane.sh}"
SMOKE_SCRIPT="${TIER2_SMOKE_SCRIPT:-$ROOT_DIR/scripts/aiida-slurm-smoke-vm.sh}"
PROFILE="${AIIDA_PROFILE:-}"
COMPUTER_LABEL="${AIIDA_COMPUTER_LABEL:-}"
WORKDIR_TEMPLATE="${AIIDA_WORKDIR_TEMPLATE:-}"
DRY_RUN=0

RUN_DIR=""
SUMMARY_JSON=""
SUMMARY_MD=""
BUNDLE_PATH=""
VALIDATE_EXIT=0
SMOKE_EXIT=3
SMOKE_SKIPPED=0
FINAL_EXIT=1
STATUS="hard-failure"
STARTED_AT_UTC=""
FINISHED_AT_UTC=""

usage() {
  cat <<'USAGE'
Usage: aiida-tier2-gate-vm.sh [options]

Run Tier-2 VM gate in one command:
  1) validate-slurm-vm-control-plane.sh --mode vm
  2) aiida-slurm-smoke-vm.sh

Options:
  --artifact-base <path>     Artifact base directory
                             (default: investigations/artifacts/gra-130)
  --run-id <id>              Deterministic run id used for artifact naming
                             (default: UTC timestamp)
  --validate-script <path>   Path to control-plane validator script
  --smoke-script <path>      Path to AiiDA->Slurm smoke script
  --profile <name>           Forwarded to smoke script
  --computer-label <label>   Forwarded to smoke script
  --workdir <path>           Forwarded to smoke script
  --dry-run                  Print planned commands without executing scripts
  -h, --help                 Show help

Exit codes:
  0  Tier-2 gate passed
  2  Blocked (runtime/prerequisite blocker)
  1  Hard failure (deterministic failure, invalid invocation, or unexpected error)
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

q() {
  printf '%q' "$1"
}

json_escape() {
  local value="$1"
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  value=${value//$'\n'/\\n}
  value=${value//$'\r'/\\r}
  value=${value//$'\t'/\\t}
  printf '%s' "$value"
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --artifact-base)
        shift
        [[ $# -gt 0 ]] || die "--artifact-base requires a value"
        ARTIFACT_BASE="$1"
        ;;
      --run-id)
        shift
        [[ $# -gt 0 ]] || die "--run-id requires a value"
        RUN_ID="$1"
        ;;
      --validate-script)
        shift
        [[ $# -gt 0 ]] || die "--validate-script requires a value"
        VALIDATE_SCRIPT="$1"
        ;;
      --smoke-script)
        shift
        [[ $# -gt 0 ]] || die "--smoke-script requires a value"
        SMOKE_SCRIPT="$1"
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

  [[ -n "$RUN_ID" ]] || die "run id must not be empty"

  ARTIFACT_BASE="$(to_abs_path "$ARTIFACT_BASE")"
  VALIDATE_SCRIPT="$(to_abs_path "$VALIDATE_SCRIPT")"
  SMOKE_SCRIPT="$(to_abs_path "$SMOKE_SCRIPT")"
}

setup_run_dir() {
  RUN_DIR="${ARTIFACT_BASE%/}/${RUN_ID}-tier2-gate-vm"
  SUMMARY_JSON="$RUN_DIR/tier2-gate-vm-summary.json"
  SUMMARY_MD="$RUN_DIR/tier2-gate-vm-summary.md"
  BUNDLE_PATH="$RUN_DIR/tier2-gate-vm-artifacts.tar.gz"

  mkdir -p "$RUN_DIR"
}

write_metadata() {
  cat >"$RUN_DIR/run-metadata.txt" <<EOF_META
run_id=$RUN_ID
artifact_base=$ARTIFACT_BASE
run_dir=$RUN_DIR
validate_script=$VALIDATE_SCRIPT
smoke_script=$SMOKE_SCRIPT
profile=${PROFILE:-<default>}
computer_label=${COMPUTER_LABEL:-<default>}
workdir_template=${WORKDIR_TEMPLATE:-<default>}
dry_run=$DRY_RUN
started_at_utc=$STARTED_AT_UTC
EOF_META
}

run_capture() {
  local log_name="$1"
  shift
  local log_path="$RUN_DIR/$log_name"

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
    log "[ok] ${log_name}"
  else
    log "[fail] ${log_name} (exit=${status})"
  fi
  return "$status"
}

collect_smoke_artifacts() {
  local smoke_root="$RUN_DIR/02-smoke-raw"
  local normalized="$RUN_DIR/02-smoke-artifacts"

  mkdir -p "$normalized"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf 'dry-run: smoke artifacts not created\n' >"$normalized/README.txt"
    return 0
  fi

  if [[ ! -d "$smoke_root" ]]; then
    printf 'smoke artifact root missing: %s\n' "$smoke_root" >"$normalized/README.txt"
    return 0
  fi

  local latest_dir
  latest_dir="$(find "$smoke_root" -mindepth 1 -maxdepth 1 -type d | LC_ALL=C sort | tail -n 1)"
  if [[ -z "$latest_dir" ]]; then
    printf 'no smoke run directory found under: %s\n' "$smoke_root" >"$normalized/README.txt"
    return 0
  fi

  cp -a "$latest_dir/." "$normalized/"
  printf 'source_smoke_run_dir=%s\n' "$latest_dir" >"$normalized/source-run.txt"
}

map_final_exit() {
  if [[ "$VALIDATE_EXIT" -eq 0 && "$SMOKE_EXIT" -eq 0 ]]; then
    FINAL_EXIT=0
    STATUS="pass"
    return 0
  fi

  if [[ "$VALIDATE_EXIT" -eq 2 || "$SMOKE_EXIT" -eq 2 ]]; then
    FINAL_EXIT=2
    STATUS="blocked"
    return 0
  fi

  FINAL_EXIT=1
  STATUS="hard-failure"
}

write_summary_json() {
  FINISHED_AT_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

  cat >"$SUMMARY_JSON" <<EOF_JSON
{
  "run_id": "$(json_escape "$RUN_ID")",
  "status": "$(json_escape "$STATUS")",
  "exit_code": $FINAL_EXIT,
  "started_at_utc": "$(json_escape "$STARTED_AT_UTC")",
  "finished_at_utc": "$(json_escape "$FINISHED_AT_UTC")",
  "validate": {
    "script": "$(json_escape "$VALIDATE_SCRIPT")",
    "exit_code": $VALIDATE_EXIT,
    "log": "01-validate-vm.log"
  },
  "smoke": {
    "script": "$(json_escape "$SMOKE_SCRIPT")",
    "exit_code": $SMOKE_EXIT,
    "skipped": $([[ "$SMOKE_SKIPPED" -eq 1 ]] && printf 'true' || printf 'false'),
    "log": "02-smoke-vm.log",
    "normalized_artifact_dir": "02-smoke-artifacts"
  },
  "artifacts": {
    "run_dir": "$(json_escape "$RUN_DIR")",
    "bundle": "tier2-gate-vm-artifacts.tar.gz",
    "summary_markdown": "tier2-gate-vm-summary.md"
  }
}
EOF_JSON
}

write_summary_md() {
  cat >"$SUMMARY_MD" <<EOF_MD
# Tier-2 Gate VM Summary

- run_id: $RUN_ID
- status: $STATUS
- final_exit: $FINAL_EXIT
- started_at_utc: $STARTED_AT_UTC
- finished_at_utc: $FINISHED_AT_UTC
- run_dir: $RUN_DIR

Checks:

- validate-slurm-vm-control-plane (--mode vm): exit $VALIDATE_EXIT
- aiida-slurm-smoke-vm: exit $SMOKE_EXIT (skipped=$SMOKE_SKIPPED)

Artifacts:

- summary json: $SUMMARY_JSON
- normalized smoke artifacts: $RUN_DIR/02-smoke-artifacts
- deterministic bundle: $BUNDLE_PATH

Exit semantics:

- 0: pass
- 2: blocked
- 1: hard failure
EOF_MD
}

create_bundle() {
  local root_name
  local parent_dir
  local tmp_bundle
  root_name="$(basename "$RUN_DIR")"
  parent_dir="$(dirname "$RUN_DIR")"
  tmp_bundle="${parent_dir}/.${root_name}.bundle.tmp.tar.gz"

  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf 'dry-run: bundle skipped\n' >"$BUNDLE_PATH"
    return 0
  fi

  if ! have_command tar; then
    die "tar is required to create artifact bundle"
  fi

  # Normalize tar metadata so bundle naming/order is deterministic.
  tar --sort=name \
    --mtime='UTC 1970-01-01' \
    --owner=0 --group=0 --numeric-owner \
    --exclude "$root_name/tier2-gate-vm-artifacts.tar.gz" \
    -C "$parent_dir" \
    -czf "$tmp_bundle" \
    "$root_name"

  mv "$tmp_bundle" "$BUNDLE_PATH"
}

run_gate() {
  local validate_cmd=(bash "$VALIDATE_SCRIPT" --mode vm)
  if run_capture "01-validate-vm.log" "${validate_cmd[@]}"; then
    VALIDATE_EXIT=0
  else
    VALIDATE_EXIT=$?
  fi

  local smoke_cmd=(bash "$SMOKE_SCRIPT" --artifact-dir "$RUN_DIR/02-smoke-raw")
  if [[ -n "$PROFILE" ]]; then
    smoke_cmd+=(--profile "$PROFILE")
  fi
  if [[ -n "$COMPUTER_LABEL" ]]; then
    smoke_cmd+=(--computer-label "$COMPUTER_LABEL")
  fi
  if [[ -n "$WORKDIR_TEMPLATE" ]]; then
    smoke_cmd+=(--workdir "$WORKDIR_TEMPLATE")
  fi

  if [[ "$VALIDATE_EXIT" -eq 0 ]]; then
    if run_capture "02-smoke-vm.log" "${smoke_cmd[@]}"; then
      SMOKE_EXIT=0
    else
      SMOKE_EXIT=$?
    fi
  else
    printf 'skip aiida smoke because validate phase failed (exit=%s)\n' "$VALIDATE_EXIT" \
      | tee "$RUN_DIR/02-smoke-vm.log"
    SMOKE_EXIT=3
    SMOKE_SKIPPED=1
  fi

  collect_smoke_artifacts
}

validate_inputs() {
  [[ -f "$VALIDATE_SCRIPT" ]] || die "validate script not found: $VALIDATE_SCRIPT"
  [[ -f "$SMOKE_SCRIPT" ]] || die "smoke script not found: $SMOKE_SCRIPT"
}

main() {
  parse_args "$@"
  validate_inputs

  STARTED_AT_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  setup_run_dir
  write_metadata

  log "Tier-2 gate VM start"
  log "  run id: $RUN_ID"
  log "  run dir: $RUN_DIR"

  run_gate
  map_final_exit
  write_summary_json
  write_summary_md
  create_bundle

  log ""
  log "Tier-2 gate VM summary"
  log "  status: $STATUS"
  log "  exit: $FINAL_EXIT"
  log "  summary json: $SUMMARY_JSON"
  log "  summary markdown: $SUMMARY_MD"
  log "  artifacts bundle: $BUNDLE_PATH"

  exit "$FINAL_EXIT"
}

main "$@"
