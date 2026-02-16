#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

MODE="docker-slurm"
ENV_FILE="$ROOT_DIR/apps/api/.env.aiida.dev"
ARTIFACT_DIR="$ROOT_DIR/investigations/artifacts/gra-84"
DIRECT_LABEL="gra84-local-direct"
SLURM_LABEL="gra84-local-slurm"
DIRECT_WORKDIR="/tmp/aiida-gra84-direct-{username}"
SLURM_WORKDIR="/tmp/aiida-gra84-slurm-{username}"
LOCAL_PROFILE_NAME="gra84-local-smoke"
LOCAL_AIIDA_PATH="$ROOT_DIR/investigations/.aiida-local-smoke/config"
LOCAL_SQLITE_PATH="$ROOT_DIR/investigations/.aiida-local-smoke/storage"
EXIT_CODE=0

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Run a minimal AiiDA smoke path for Docker/Slurm validation.

Modes:
  docker-slurm   Bootstraps core.psql_dos profile via docker, runs direct smoke,
                 then attempts slurm smoke. If slurm is unavailable, exits 2.
  local-only     Creates a local core.sqlite_dos profile and runs direct smoke.

Options:
  --mode <mode>           One of: docker-slurm, local-only (default: docker-slurm)
  --env-file <path>       Env file for docker-slurm mode (default: apps/api/.env.aiida.dev)
  --artifact-dir <path>   Directory for logs (default: investigations/artifacts/gra-84)
  -h, --help              Show help
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
      --mode)
        shift
        [[ $# -gt 0 ]] || die "--mode requires a value"
        MODE="$1"
        ;;
      --env-file)
        shift
        [[ $# -gt 0 ]] || die "--env-file requires a value"
        ENV_FILE="$1"
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

  if [[ "$MODE" != "docker-slurm" && "$MODE" != "local-only" ]]; then
    die "unsupported mode: $MODE"
  fi

  ENV_FILE="$(to_abs_path "$ENV_FILE")"
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
  local name="$1"
  local log_path="$RUN_DIR/${name}.log"
  if grep -Eq '\\[Failed\\]|tests failed' "$log_path"; then
    return 1
  fi
  return 0
}

profile_exists() {
  local aiida_path="$1"
  local profile_name="$2"
  AIIDA_PATH="$aiida_path" uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi profile show "$profile_name" >/dev/null 2>&1
}

computer_exists() {
  local aiida_path="$1"
  local profile_name="$2"
  local label="$3"
  AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer show "$label" >/dev/null 2>&1
}

ensure_direct_computer() {
  local aiida_path="$1"
  local profile_name="$2"

  if computer_exists "$aiida_path" "$profile_name" "$DIRECT_LABEL"; then
    log "direct computer already exists: ${DIRECT_LABEL}"
    return 0
  fi

  AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer setup --non-interactive \
    --label "$DIRECT_LABEL" \
    --hostname localhost \
    --transport core.local \
    --scheduler core.direct \
    --work-dir "$DIRECT_WORKDIR"

  AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer configure core.local --non-interactive \
    --safe-interval 0 "$DIRECT_LABEL"
}

ensure_slurm_computer() {
  local aiida_path="$1"
  local profile_name="$2"

  if computer_exists "$aiida_path" "$profile_name" "$SLURM_LABEL"; then
    log "slurm computer already exists: ${SLURM_LABEL}"
    return 0
  fi

  AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer setup --non-interactive \
    --label "$SLURM_LABEL" \
    --hostname localhost \
    --transport core.local \
    --scheduler core.slurm \
    --work-dir "$SLURM_WORKDIR"

  AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer configure core.local --non-interactive \
    --safe-interval 0 "$SLURM_LABEL"
}

run_direct_smoke() {
  local aiida_path="$1"
  local profile_name="$2"

  ensure_direct_computer "$aiida_path" "$profile_name"
  run_capture "aiida-computer-test-direct" \
    env AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer test "$DIRECT_LABEL"
  computer_test_clean "aiida-computer-test-direct"
}

run_slurm_smoke() {
  local aiida_path="$1"
  local profile_name="$2"

  ensure_slurm_computer "$aiida_path" "$profile_name"

  if ! run_capture "slurm-sinfo" sinfo; then
    log "blocker: slurm control plane not reachable (sinfo failed)."
  fi

  if run_capture "aiida-computer-test-slurm" \
    env AIIDA_PATH="$aiida_path" AIIDA_PROFILE="$profile_name" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi computer test "$SLURM_LABEL" && \
    computer_test_clean "aiida-computer-test-slurm"; then
    return 0
  fi

  log "blocker: AiiDA slurm computer smoke failed."
  log "fallback: run local scheduler smoke mode:"
  log "  scripts/aiida-slurm-smoke-docker.sh --mode local-only"
  return 2
}

load_env_file() {
  [[ -f "$ENV_FILE" ]] || die "missing env file: $ENV_FILE"
  set -a
  # shellcheck disable=SC1090
  source "$ENV_FILE"
  set +a
}

prepare_docker_profile() {
  have_command docker || die "docker is required in docker-slurm mode"

  load_env_file

  local aiida_config_dir="${AIIDA_CONFIG_DIR:-.just-runtime/aiida/config}"
  local aiida_profile="${AIIDA_PROFILE:-chem-model-dev}"
  local aiida_repository_dir="${AIIDA_REPOSITORY_DIR:-.just-runtime/aiida/repository}"
  local aiida_pghost="${AIIDA_PGHOST:-127.0.0.1}"
  local aiida_pgport="${AIIDA_PGPORT:-5432}"
  local aiida_pgdatabase="${AIIDA_PGDATABASE:-aiida_db}"
  local aiida_pguser="${AIIDA_PGUSER:-aiida}"
  local aiida_pgpassword="${AIIDA_PGPASSWORD:-aiida-dev-only}"
  local aiida_default_email="${AIIDA_DEFAULT_USER_EMAIL:-dev@chem-model-edit.local}"
  local aiida_default_first_name="${AIIDA_DEFAULT_USER_FIRST_NAME:-Chem}"
  local aiida_default_last_name="${AIIDA_DEFAULT_USER_LAST_NAME:-Model}"
  local aiida_default_institution="${AIIDA_DEFAULT_USER_INSTITUTION:-Chem Model Edit}"
  local aiida_path_abs
  local aiida_repository_abs
  local aiida_repository_uri
  aiida_path_abs="$(to_abs_path "$aiida_config_dir")"
  aiida_repository_abs="$(to_abs_path "$aiida_repository_dir")"
  aiida_repository_uri="file://$aiida_repository_abs"

  run_capture "bootstrap-aiida-runtime" \
    "$ROOT_DIR/scripts/bootstrap-aiida-runtime.sh" --apply --env-file "$ENV_FILE"

  if profile_exists "$aiida_path_abs" "$aiida_profile"; then
    log "docker profile already exists: ${aiida_profile}"
  else
    run_capture "setup-docker-psql-profile" \
      env AIIDA_PATH="$aiida_path_abs" \
      uv run --project "$ROOT_DIR/apps/api" --group aiida \
      verdi profile setup core.psql_dos --non-interactive \
      --profile-name "$aiida_profile" \
      --email "$aiida_default_email" \
      --first-name "$aiida_default_first_name" \
      --last-name "$aiida_default_last_name" \
      --institution "$aiida_default_institution" \
      --use-rabbitmq \
      --database-hostname "$aiida_pghost" \
      --database-port "$aiida_pgport" \
      --database-name "$aiida_pgdatabase" \
      --database-username "$aiida_pguser" \
      --database-password "$aiida_pgpassword" \
      --repository-uri "$aiida_repository_uri"
  fi

  printf '%s\n' "$aiida_path_abs" >"$RUN_DIR/docker_aiida_path.txt"
  printf '%s\n' "$aiida_profile" >"$RUN_DIR/docker_profile_name.txt"
}

prepare_local_profile() {
  mkdir -p "$LOCAL_AIIDA_PATH" "$LOCAL_SQLITE_PATH"

  if profile_exists "$LOCAL_AIIDA_PATH" "$LOCAL_PROFILE_NAME"; then
    log "local sqlite profile already exists: ${LOCAL_PROFILE_NAME}"
    return 0
  fi

  run_capture "setup-local-sqlite-profile" \
    env AIIDA_PATH="$LOCAL_AIIDA_PATH" \
    uv run --project "$ROOT_DIR/apps/api" --group aiida \
    verdi profile setup core.sqlite_dos --non-interactive \
    --profile-name "$LOCAL_PROFILE_NAME" \
    --email dev@chem-model-edit.local \
    --first-name Chem \
    --last-name Model \
    --institution "Chem Model Edit" \
    --no-use-rabbitmq \
    --filepath "$LOCAL_SQLITE_PATH"
}

main() {
  parse_args "$@"

  have_command uv || die "uv is required"
  mkdir -p "$ARTIFACT_DIR"
  RUN_DIR="$ARTIFACT_DIR/$(date +%Y%m%d-%H%M%S)-${MODE}"
  mkdir -p "$RUN_DIR"

  log "AiiDA smoke start"
  log "  mode: ${MODE}"
  log "  run dir: ${RUN_DIR}"

  if [[ "$MODE" == "docker-slurm" ]]; then
    prepare_docker_profile
    local aiida_path
    local profile_name
    aiida_path="$(cat "$RUN_DIR/docker_aiida_path.txt")"
    profile_name="$(cat "$RUN_DIR/docker_profile_name.txt")"

    run_direct_smoke "$aiida_path" "$profile_name" || EXIT_CODE=1
    if [[ "$EXIT_CODE" -eq 0 ]]; then
      run_slurm_smoke "$aiida_path" "$profile_name" || EXIT_CODE=$?
    fi
  else
    prepare_local_profile
    run_direct_smoke "$LOCAL_AIIDA_PATH" "$LOCAL_PROFILE_NAME" || EXIT_CODE=1
  fi

  log ""
  log "Smoke summary"
  log "  mode: ${MODE}"
  log "  artifacts: ${RUN_DIR}"

  if [[ "$EXIT_CODE" -eq 0 ]]; then
    log "  result: success"
  elif [[ "$EXIT_CODE" -eq 2 ]]; then
    log "  result: blocked at slurm stage (fallback available)"
  else
    log "  result: failed"
  fi

  exit "$EXIT_CODE"
}

main "$@"
