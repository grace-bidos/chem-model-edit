#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ARTIFACT_DIR="$ROOT_DIR/investigations/artifacts/gra-105"
RUN_LABEL="local"
EXIT_CODE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Run a lightweight local/dev smoke harness for submit -> event -> projection lifecycle.

Options:
  --artifact-dir <path>   Directory for smoke artifacts
                          (default: investigations/artifacts/gra-105)
  --run-label <label>     Label added to run directory name (default: local)
  -h, --help              Show help

Notes:
  - This harness is Slurm-independent and does not invoke Slurm commands.
  - It runs only targeted API contract/lifecycle tests.
USAGE
}

log() {
  printf '%s\n' "$*"
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
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
      --artifact-dir)
        shift
        [[ $# -gt 0 ]] || die "--artifact-dir requires a value"
        ARTIFACT_DIR="$1"
        ;;
      --run-label)
        shift
        [[ $# -gt 0 ]] || die "--run-label requires a value"
        RUN_LABEL="$1"
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
  RUN_LABEL="$(printf '%s' "$RUN_LABEL" | tr -cs 'a-zA-Z0-9_.-' '-')"
}

require_command() {
  local name="$1"
  command -v "$name" >/dev/null 2>&1 || die "required command not found: $name"
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

write_manifest() {
  cat > "$RUN_DIR/manifest.txt" <<MANIFEST
smoke: submit->event->projection
started_at_utc: $STARTED_AT_UTC
run_dir: $RUN_DIR
host: ${HOSTNAME:-unknown}
branch: $(git -C "$ROOT_DIR" rev-parse --abbrev-ref HEAD)
commit: $(git -C "$ROOT_DIR" rev-parse HEAD)
MANIFEST
}

main() {
  parse_args "$@"
  require_command uv
  require_command git

  STARTED_AT_UTC="$(date -u +"%Y%m%dT%H%M%SZ")"
  RUN_DIR="$ARTIFACT_DIR/${STARTED_AT_UTC}-${RUN_LABEL}"
  mkdir -p "$RUN_DIR"

  write_manifest

  log "Starting lifecycle smoke: submit -> event -> projection"
  log "Artifacts: $RUN_DIR"
  log "This smoke harness is Slurm-independent."

  if ! run_capture "pytest-submit-event-projection" \
    uv run --project "$ROOT_DIR/apps/api" pytest \
      "$ROOT_DIR/apps/api/tests/test_zpe_compute_results.py::test_submit_result_idempotent" \
      "$ROOT_DIR/apps/api/tests/test_zpe_compute_results.py::test_submit_result_event_id_replay_is_idempotent_and_payload_safe" \
      "$ROOT_DIR/apps/api/tests/test_zpe_compute_results.py::test_submit_failure_event_id_replay_is_idempotent_and_payload_safe" \
      "$ROOT_DIR/apps/api/tests/test_zpe_submit_event_lifecycle_contract.py::test_result_request_normalizes_aiida_runtime_execution_event_shape" \
      "$ROOT_DIR/apps/api/tests/test_zpe_submit_event_lifecycle_contract.py::test_failed_request_normalizes_aiida_runtime_execution_event_shape" \
      "$ROOT_DIR/apps/api/tests/test_convex_event_relay.py::test_projection_as_dict_matches_contract" \
      "$ROOT_DIR/apps/api/tests/test_convex_event_relay.py::test_http_dispatcher_posts_projection_with_idempotency_header"; then
    EXIT_CODE=1
  fi

  if [[ "$EXIT_CODE" -eq 0 ]]; then
    log ""
    log "Smoke result: PASS"
  else
    log ""
    log "Smoke result: FAIL"
  fi

  log "Logs and artifacts: $RUN_DIR"
  exit "$EXIT_CODE"
}

main "$@"
