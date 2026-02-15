#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_ENV_FILE="$ROOT_DIR/apps/api/.env.aiida.dev"
EXAMPLE_ENV_FILE="$ROOT_DIR/apps/api/.env.aiida.dev.example"

DRY_RUN=1
COPY_ENV=0
INIT_PROFILE=0
ENV_FILE="$DEFAULT_ENV_FILE"

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Bootstrap a minimal AiiDA runtime (PostgreSQL + RabbitMQ) for dev/staging.
Safe by default: runs in dry-run mode unless --apply is set.

Options:
  --apply                 Execute commands (default: dry-run only)
  --copy-env              Copy apps/api/.env.aiida.dev.example to env file if missing
  --env-file <path>       Environment file path (default: apps/api/.env.aiida.dev)
  --init-profile          Also run 'verdi profile setup core.psql_dos' if profile is missing
  -h, --help              Show this help
EOF
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

run_cmd() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run]'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

to_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

docker_network_exists() {
  docker network inspect "$1" >/dev/null 2>&1
}

docker_container_exists() {
  docker ps -a --format '{{.Names}}' | grep -Fxq "$1"
}

docker_container_running() {
  docker ps --format '{{.Names}}' | grep -Fxq "$1"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)
      DRY_RUN=0
      ;;
    --copy-env)
      COPY_ENV=1
      ;;
    --env-file)
      shift
      [[ $# -gt 0 ]] || die "--env-file requires a value"
      ENV_FILE="$1"
      ;;
    --init-profile)
      INIT_PROFILE=1
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

if [[ "$ENV_FILE" != /* ]]; then
  ENV_FILE="$ROOT_DIR/$ENV_FILE"
fi

if [[ ! -f "$ENV_FILE" && "$COPY_ENV" -eq 1 ]]; then
  [[ -f "$EXAMPLE_ENV_FILE" ]] || die "missing example env file: $EXAMPLE_ENV_FILE"
  run_cmd cp "$EXAMPLE_ENV_FILE" "$ENV_FILE"
fi

[[ -f "$ENV_FILE" ]] || die "missing env file: $ENV_FILE (tip: pass --copy-env)"

set -a
source "$ENV_FILE"
set +a

AIIDA_PROFILE="${AIIDA_PROFILE:-chem-model-dev}"
AIIDA_ENV_NAME="${AIIDA_ENV_NAME:-dev}"
AIIDA_CONFIG_DIR="${AIIDA_CONFIG_DIR:-.just-runtime/aiida/config}"
AIIDA_REPOSITORY_DIR="${AIIDA_REPOSITORY_DIR:-.just-runtime/aiida/repository}"
AIIDA_DOCKER_NETWORK="${AIIDA_DOCKER_NETWORK:-chem-model-aiida-net}"
AIIDA_PG_CONTAINER_NAME="${AIIDA_PG_CONTAINER_NAME:-chem-model-aiida-pg}"
AIIDA_RMQ_CONTAINER_NAME="${AIIDA_RMQ_CONTAINER_NAME:-chem-model-aiida-rmq}"
AIIDA_PGHOST="${AIIDA_PGHOST:-127.0.0.1}"
AIIDA_PGPORT="${AIIDA_PGPORT:-5432}"
AIIDA_PGDATABASE="${AIIDA_PGDATABASE:-aiida_db}"
AIIDA_PGUSER="${AIIDA_PGUSER:-aiida}"
AIIDA_PGPASSWORD="${AIIDA_PGPASSWORD:-aiida-dev-only}"
AIIDA_BROKER_PROTOCOL="${AIIDA_BROKER_PROTOCOL:-amqp}"
AIIDA_BROKER_HOST="${AIIDA_BROKER_HOST:-127.0.0.1}"
AIIDA_BROKER_PORT="${AIIDA_BROKER_PORT:-5672}"
AIIDA_BROKER_USER="${AIIDA_BROKER_USER:-guest}"
AIIDA_BROKER_PASSWORD="${AIIDA_BROKER_PASSWORD:-guest}"
AIIDA_RMQ_MANAGEMENT_PORT="${AIIDA_RMQ_MANAGEMENT_PORT:-15672}"
AIIDA_DEFAULT_USER_EMAIL="${AIIDA_DEFAULT_USER_EMAIL:-dev@chem-model-edit.local}"
AIIDA_DEFAULT_USER_FIRST_NAME="${AIIDA_DEFAULT_USER_FIRST_NAME:-Chem}"
AIIDA_DEFAULT_USER_LAST_NAME="${AIIDA_DEFAULT_USER_LAST_NAME:-Model}"
AIIDA_DEFAULT_USER_INSTITUTION="${AIIDA_DEFAULT_USER_INSTITUTION:-Chem Model Edit}"

AIIDA_CONFIG_DIR_ABS="$(to_abs_path "$AIIDA_CONFIG_DIR")"
AIIDA_REPOSITORY_DIR_ABS="$(to_abs_path "$AIIDA_REPOSITORY_DIR")"
AIIDA_RUNTIME_DIR_ABS="$(dirname "$AIIDA_CONFIG_DIR_ABS")"
AIIDA_PGDATA_DIR_ABS="$AIIDA_RUNTIME_DIR_ABS/postgres-data"
AIIDA_RMQDATA_DIR_ABS="$AIIDA_RUNTIME_DIR_ABS/rabbitmq-data"
AIIDA_REPOSITORY_URI="${AIIDA_REPOSITORY_URI:-file://$AIIDA_REPOSITORY_DIR_ABS}"

if [[ "$DRY_RUN" -eq 0 ]]; then
  have_command docker || die "docker is required with --apply"
  if [[ "$INIT_PROFILE" -eq 1 ]]; then
    have_command uv || die "uv is required with --init-profile"
  fi
fi

log "AiiDA runtime bootstrap"
log "  Mode:              ${AIIDA_ENV_NAME}"
log "  Dry run:           $([[ "$DRY_RUN" -eq 1 ]] && echo yes || echo no)"
log "  Env file:          ${ENV_FILE}"
log "  AiiDA profile:     ${AIIDA_PROFILE}"
log "  AiiDA config dir:  ${AIIDA_CONFIG_DIR_ABS}"
log "  AiiDA repo dir:    ${AIIDA_REPOSITORY_DIR_ABS}"
log "  Docker network:    ${AIIDA_DOCKER_NETWORK}"
log "  Postgres:          ${AIIDA_PG_CONTAINER_NAME} (${AIIDA_PGHOST}:${AIIDA_PGPORT}/${AIIDA_PGDATABASE})"
log "  RabbitMQ:          ${AIIDA_RMQ_CONTAINER_NAME} (${AIIDA_BROKER_HOST}:${AIIDA_BROKER_PORT})"

run_cmd mkdir -p "$AIIDA_RUNTIME_DIR_ABS" "$AIIDA_CONFIG_DIR_ABS" "$AIIDA_REPOSITORY_DIR_ABS" \
  "$AIIDA_PGDATA_DIR_ABS" "$AIIDA_RMQDATA_DIR_ABS"

if have_command docker; then
  if docker_network_exists "$AIIDA_DOCKER_NETWORK"; then
    log "docker network exists: ${AIIDA_DOCKER_NETWORK}"
  else
    run_cmd docker network create "$AIIDA_DOCKER_NETWORK"
  fi

  if docker_container_exists "$AIIDA_PG_CONTAINER_NAME"; then
    if docker_container_running "$AIIDA_PG_CONTAINER_NAME"; then
      log "postgres container already running: ${AIIDA_PG_CONTAINER_NAME}"
    else
      run_cmd docker start "$AIIDA_PG_CONTAINER_NAME"
    fi
  else
    run_cmd docker run -d --name "$AIIDA_PG_CONTAINER_NAME" \
      --network "$AIIDA_DOCKER_NETWORK" \
      -e POSTGRES_DB="$AIIDA_PGDATABASE" \
      -e POSTGRES_USER="$AIIDA_PGUSER" \
      -e POSTGRES_PASSWORD="$AIIDA_PGPASSWORD" \
      -p "${AIIDA_PGPORT}:5432" \
      -v "${AIIDA_PGDATA_DIR_ABS}:/var/lib/postgresql/data" \
      postgres:16
  fi

  if docker_container_exists "$AIIDA_RMQ_CONTAINER_NAME"; then
    if docker_container_running "$AIIDA_RMQ_CONTAINER_NAME"; then
      log "rabbitmq container already running: ${AIIDA_RMQ_CONTAINER_NAME}"
    else
      run_cmd docker start "$AIIDA_RMQ_CONTAINER_NAME"
    fi
  else
    run_cmd docker run -d --name "$AIIDA_RMQ_CONTAINER_NAME" \
      --network "$AIIDA_DOCKER_NETWORK" \
      -p "${AIIDA_BROKER_PORT}:5672" \
      -p "${AIIDA_RMQ_MANAGEMENT_PORT}:15672" \
      -v "${AIIDA_RMQDATA_DIR_ABS}:/var/lib/rabbitmq" \
      rabbitmq:3.13-management
  fi
else
  log "warning: docker not found on PATH; container checks were skipped."
fi

if [[ "$INIT_PROFILE" -eq 1 ]]; then
  profile_exists=1
  if [[ "$DRY_RUN" -eq 0 ]]; then
    if AIIDA_PATH="$AIIDA_CONFIG_DIR_ABS" uv run --project "$ROOT_DIR/apps/api" --group aiida \
      verdi profile show "$AIIDA_PROFILE" >/dev/null 2>&1; then
      profile_exists=0
    fi
  fi

  if [[ "$DRY_RUN" -eq 1 || "$profile_exists" -ne 0 ]]; then
    run_cmd env "AIIDA_PATH=$AIIDA_CONFIG_DIR_ABS" \
      uv run --project "$ROOT_DIR/apps/api" --group aiida \
      verdi profile setup core.psql_dos \
      --non-interactive \
      --profile "$AIIDA_PROFILE" \
      --email "$AIIDA_DEFAULT_USER_EMAIL" \
      --first-name "$AIIDA_DEFAULT_USER_FIRST_NAME" \
      --last-name "$AIIDA_DEFAULT_USER_LAST_NAME" \
      --institution "$AIIDA_DEFAULT_USER_INSTITUTION" \
      --db-host "$AIIDA_PGHOST" \
      --db-port "$AIIDA_PGPORT" \
      --db-name "$AIIDA_PGDATABASE" \
      --db-username "$AIIDA_PGUSER" \
      --db-password "$AIIDA_PGPASSWORD" \
      --broker-protocol "$AIIDA_BROKER_PROTOCOL" \
      --broker-host "$AIIDA_BROKER_HOST" \
      --broker-port "$AIIDA_BROKER_PORT" \
      --broker-username "$AIIDA_BROKER_USER" \
      --broker-password "$AIIDA_BROKER_PASSWORD" \
      --repository "$AIIDA_REPOSITORY_URI"
  else
    log "AiiDA profile already exists: ${AIIDA_PROFILE}"
  fi
fi

log ""
log "Next steps:"
log "  1) Install dependencies:"
log "     uv sync --project \"$ROOT_DIR/apps/api\" --group aiida"
log "  2) Optional profile setup (one time):"
log "     $0 --apply --env-file \"$ENV_FILE\" --init-profile"
log "  3) Verify runtime:"
log "     AIIDA_PATH=\"$AIIDA_CONFIG_DIR_ABS\" uv run --project \"$ROOT_DIR/apps/api\" --group aiida verdi status"
