#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
DEFAULT_ENV_FILE="$ROOT_DIR/ops/aiida-vm/aiida-vm.env"
EXAMPLE_ENV_FILE="$ROOT_DIR/ops/aiida-vm/aiida-vm.env.example"

DRY_RUN=1
COPY_ENV=0
INIT_PROFILE=0
SANITY_CHECK=0
ENV_FILE="$DEFAULT_ENV_FILE"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Bootstrap AiiDA runtime on a VM (PostgreSQL + RabbitMQ + AiiDA profile).
Safe by default: runs in dry-run mode unless --apply is set.

Options:
  --apply                 Execute commands (default: dry-run only)
  --copy-env              Copy ops/aiida-vm/aiida-vm.env.example if env file is missing
  --env-file <path>       Environment file path (default: ops/aiida-vm/aiida-vm.env)
  --init-profile          Create AiiDA profile if missing
  --sanity-check          Run "verdi profile list" and "verdi status"
  -h, --help              Show this help
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

run_cmd() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run]'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  "$@"
}

run_cmd_as_postgres() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run] sudo -u postgres'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  sudo -u postgres "$@"
}

run_cmd_with_sudo() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    printf '[dry-run] sudo'
    printf ' %q' "$@"
    printf '\n'
    return 0
  fi
  sudo "$@"
}

to_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

wait_for_tcp() {
  local host="$1"
  local port="$2"
  local timeout_seconds="$3"
  local start_ts now_ts
  start_ts="$(date +%s)"
  while true; do
    if python3 - "$host" "$port" <<'PY' >/dev/null 2>&1
import socket
import sys

host = sys.argv[1]
port = int(sys.argv[2])

sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
sock.settimeout(1.0)
try:
    sock.connect((host, port))
except OSError:
    raise SystemExit(1)
finally:
    sock.close()
PY
    then
      return 0
    fi
    now_ts="$(date +%s)"
    if (( now_ts - start_ts >= timeout_seconds )); then
      return 1
    fi
    sleep 1
  done
}

sql_quote() {
  local value="$1"
  value="$(printf '%s' "$value" | sed "s/'/''/g")"
  printf "'%s'" "$value"
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
    --sanity-check)
      SANITY_CHECK=1
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

if [[ ! -f "$ENV_FILE" && "$DRY_RUN" -eq 1 && "$COPY_ENV" -eq 1 ]]; then
  log "env file is missing in dry-run; using example template for preview: $EXAMPLE_ENV_FILE"
  ENV_FILE="$EXAMPLE_ENV_FILE"
fi

[[ -f "$ENV_FILE" ]] || die "missing env file: $ENV_FILE (tip: pass --copy-env)"

set -a
source "$ENV_FILE"
set +a

AIIDA_VM_APT_UPDATE="${AIIDA_VM_APT_UPDATE:-1}"
AIIDA_VM_INSTALL_PACKAGES="${AIIDA_VM_INSTALL_PACKAGES:-1}"
AIIDA_VM_ENABLE_SERVICES="${AIIDA_VM_ENABLE_SERVICES:-1}"
AIIDA_VM_WAIT_TIMEOUT_SECONDS="${AIIDA_VM_WAIT_TIMEOUT_SECONDS:-60}"

AIIDA_PGHOST="${AIIDA_PGHOST:-127.0.0.1}"
AIIDA_PGPORT="${AIIDA_PGPORT:-5432}"
AIIDA_PGDATABASE="${AIIDA_PGDATABASE:-aiida_db}"
AIIDA_PGUSER="${AIIDA_PGUSER:-aiida}"
AIIDA_PGPASSWORD="${AIIDA_PGPASSWORD:-aiida-dev-only}"

AIIDA_BROKER_PROTOCOL="${AIIDA_BROKER_PROTOCOL:-amqp}"
AIIDA_BROKER_HOST="${AIIDA_BROKER_HOST:-127.0.0.1}"
AIIDA_BROKER_PORT="${AIIDA_BROKER_PORT:-5672}"
AIIDA_BROKER_USER="${AIIDA_BROKER_USER:-aiida}"
AIIDA_BROKER_PASSWORD="${AIIDA_BROKER_PASSWORD:-aiida-dev-only}"
AIIDA_BROKER_VHOST="${AIIDA_BROKER_VHOST:-aiida}"

AIIDA_PROFILE="${AIIDA_PROFILE:-chem-model-vm}"
AIIDA_CONFIG_DIR="${AIIDA_CONFIG_DIR:-.just-runtime/aiida-vm/config}"
AIIDA_REPOSITORY_DIR="${AIIDA_REPOSITORY_DIR:-.just-runtime/aiida-vm/repository}"
AIIDA_DEFAULT_USER_EMAIL="${AIIDA_DEFAULT_USER_EMAIL:-dev@chem-model-edit.local}"
AIIDA_DEFAULT_USER_FIRST_NAME="${AIIDA_DEFAULT_USER_FIRST_NAME:-Chem}"
AIIDA_DEFAULT_USER_LAST_NAME="${AIIDA_DEFAULT_USER_LAST_NAME:-Model}"
AIIDA_DEFAULT_USER_INSTITUTION="${AIIDA_DEFAULT_USER_INSTITUTION:-Chem Model Edit}"
AIIDA_UV_PROJECT="${AIIDA_UV_PROJECT:-apps/api}"

AIIDA_CONFIG_DIR_ABS="$(to_abs_path "$AIIDA_CONFIG_DIR")"
AIIDA_REPOSITORY_DIR_ABS="$(to_abs_path "$AIIDA_REPOSITORY_DIR")"
AIIDA_REPOSITORY_URI="file://$AIIDA_REPOSITORY_DIR_ABS"
AIIDA_UV_PROJECT_ABS="$(to_abs_path "$AIIDA_UV_PROJECT")"

if [[ "$DRY_RUN" -eq 0 ]]; then
  have_command sudo || die "sudo is required with --apply"
  have_command uv || die "uv is required with --apply"
  have_command python3 || die "python3 is required with --apply"
fi

log "AiiDA VM bootstrap"
log "  Dry run:           $([[ "$DRY_RUN" -eq 1 ]] && echo yes || echo no)"
log "  Env file:          $ENV_FILE"
log "  AiiDA profile:     $AIIDA_PROFILE"
log "  AiiDA config dir:  $AIIDA_CONFIG_DIR_ABS"
log "  AiiDA repo dir:    $AIIDA_REPOSITORY_DIR_ABS"
log "  PostgreSQL:        $AIIDA_PGHOST:$AIIDA_PGPORT/$AIIDA_PGDATABASE"
log "  RabbitMQ:          $AIIDA_BROKER_HOST:$AIIDA_BROKER_PORT vhost=$AIIDA_BROKER_VHOST"

run_cmd mkdir -p "$AIIDA_CONFIG_DIR_ABS" "$AIIDA_REPOSITORY_DIR_ABS"

if [[ "$AIIDA_VM_APT_UPDATE" = "1" ]]; then
  run_cmd_with_sudo apt-get update
fi

if [[ "$AIIDA_VM_INSTALL_PACKAGES" = "1" ]]; then
  run_cmd_with_sudo apt-get install -y postgresql rabbitmq-server
fi

if [[ "$AIIDA_VM_ENABLE_SERVICES" = "1" ]]; then
  run_cmd_with_sudo systemctl enable --now postgresql rabbitmq-server
fi

if [[ "$DRY_RUN" -eq 0 ]]; then
  log "waiting for postgres readiness: $AIIDA_PGHOST:$AIIDA_PGPORT"
  wait_for_tcp "$AIIDA_PGHOST" "$AIIDA_PGPORT" "$AIIDA_VM_WAIT_TIMEOUT_SECONDS" \
    || die "postgres did not become ready within timeout"
  log "waiting for rabbitmq readiness: $AIIDA_BROKER_HOST:$AIIDA_BROKER_PORT"
  wait_for_tcp "$AIIDA_BROKER_HOST" "$AIIDA_BROKER_PORT" "$AIIDA_VM_WAIT_TIMEOUT_SECONDS" \
    || die "rabbitmq did not become ready within timeout"
fi

DB_USER_SQL="DO \\$\\$ BEGIN\\n  IF NOT EXISTS (SELECT 1 FROM pg_roles WHERE rolname = $(sql_quote "$AIIDA_PGUSER")) THEN\\n    CREATE USER \\\"$AIIDA_PGUSER\\\" WITH PASSWORD $(sql_quote "$AIIDA_PGPASSWORD");\\n  ELSE\\n    ALTER USER \\\"$AIIDA_PGUSER\\\" WITH PASSWORD $(sql_quote "$AIIDA_PGPASSWORD");\\n  END IF;\\nEND \\$\\$;"

run_cmd_as_postgres psql -v ON_ERROR_STOP=1 -d postgres -c "$DB_USER_SQL"
if [[ "$DRY_RUN" -eq 1 ]]; then
  run_cmd_as_postgres psql -v ON_ERROR_STOP=1 -d postgres -tAc \
    "SELECT 1 FROM pg_database WHERE datname = $(sql_quote "$AIIDA_PGDATABASE")"
  run_cmd_as_postgres createdb -O "$AIIDA_PGUSER" "$AIIDA_PGDATABASE"
else
  if ! sudo -u postgres psql -d postgres -tAc \
    "SELECT 1 FROM pg_database WHERE datname = $(sql_quote "$AIIDA_PGDATABASE")" | grep -q 1; then
    sudo -u postgres createdb -O "$AIIDA_PGUSER" "$AIIDA_PGDATABASE"
  fi
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  run_cmd_with_sudo rabbitmqctl list_vhosts
  run_cmd_with_sudo rabbitmqctl add_vhost "$AIIDA_BROKER_VHOST"
  run_cmd_with_sudo rabbitmqctl list_users
  run_cmd_with_sudo rabbitmqctl add_user "$AIIDA_BROKER_USER" "$AIIDA_BROKER_PASSWORD"
  run_cmd_with_sudo rabbitmqctl set_permissions -p "$AIIDA_BROKER_VHOST" "$AIIDA_BROKER_USER" ".*" ".*" ".*"
else
  if ! sudo rabbitmqctl list_vhosts -q | grep -Fxq "$AIIDA_BROKER_VHOST"; then
    sudo rabbitmqctl add_vhost "$AIIDA_BROKER_VHOST"
  fi
  if ! sudo rabbitmqctl list_users -q | awk '{print $1}' | grep -Fxq "$AIIDA_BROKER_USER"; then
    sudo rabbitmqctl add_user "$AIIDA_BROKER_USER" "$AIIDA_BROKER_PASSWORD"
  else
    sudo rabbitmqctl change_password "$AIIDA_BROKER_USER" "$AIIDA_BROKER_PASSWORD"
  fi
  sudo rabbitmqctl set_permissions -p "$AIIDA_BROKER_VHOST" "$AIIDA_BROKER_USER" ".*" ".*" ".*"
fi

if [[ "$INIT_PROFILE" -eq 1 ]]; then
  profile_exists=1
  if [[ "$DRY_RUN" -eq 0 ]]; then
    if AIIDA_PATH="$AIIDA_CONFIG_DIR_ABS" uv run --project "$AIIDA_UV_PROJECT_ABS" --group aiida \
      verdi profile show "$AIIDA_PROFILE" >/dev/null 2>&1; then
      profile_exists=0
    fi
  fi

  if [[ "$DRY_RUN" -eq 1 || "$profile_exists" -ne 0 ]]; then
    run_cmd env "AIIDA_PATH=$AIIDA_CONFIG_DIR_ABS" \
      uv run --project "$AIIDA_UV_PROJECT_ABS" --group aiida \
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
      --broker-virtual-host "$AIIDA_BROKER_VHOST" \
      --repository "$AIIDA_REPOSITORY_URI"
  else
    log "AiiDA profile already exists: $AIIDA_PROFILE"
  fi
fi

if [[ "$SANITY_CHECK" -eq 1 ]]; then
  run_cmd env "AIIDA_PATH=$AIIDA_CONFIG_DIR_ABS" \
    uv run --project "$AIIDA_UV_PROJECT_ABS" --group aiida \
    verdi profile list
  run_cmd env "AIIDA_PATH=$AIIDA_CONFIG_DIR_ABS" \
    uv run --project "$AIIDA_UV_PROJECT_ABS" --group aiida \
    verdi status
fi

log ""
log "Done. Suggested flow:"
log "  1) Dry-run"
log "     $0 --copy-env --env-file ops/aiida-vm/aiida-vm.env --init-profile --sanity-check"
log "  2) Apply"
log "     $0 --apply --env-file ops/aiida-vm/aiida-vm.env --init-profile --sanity-check"
