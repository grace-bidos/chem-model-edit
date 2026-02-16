#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${SCHEMATHESIS_PYTHON_BIN:-}"
if [ -z "${PYTHON_BIN}" ]; then
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_BIN="python3"
  elif command -v python >/dev/null 2>&1; then
    PYTHON_BIN="python"
  else
    echo "Neither python3 nor python is available in PATH." >&2
    echo "Set SCHEMATHESIS_PYTHON_BIN to a valid interpreter path." >&2
    exit 1
  fi
fi

HOST="${SCHEMATHESIS_HOST:-127.0.0.1}"
if [ -n "${SCHEMATHESIS_PORT:-}" ]; then
  PORT="${SCHEMATHESIS_PORT}"
else
  PORT="$(
    "${PYTHON_BIN}" - <<'PY'
import socket

with socket.socket() as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
  )"
fi
BASE_URL="${SCHEMATHESIS_BASE_URL:-http://${HOST}:${PORT}}"
SCHEMA_URL="${SCHEMATHESIS_SCHEMA_URL:-${BASE_URL}/api/openapi.json}"
SCHEMATHESIS_MODE="${SCHEMATHESIS_MODE:-smoke}"
SCHEMATHESIS_SEED="${SCHEMATHESIS_SEED:-137}"
SCHEMATHESIS_STARTUP_TIMEOUT="${SCHEMATHESIS_STARTUP_TIMEOUT:-30}"

SERVER_LOG="$(mktemp -t schemathesis-api.XXXXXX.log)"
uvicorn main:app --host "${HOST}" --port "${PORT}" --log-level warning >"${SERVER_LOG}" 2>&1 &
server_pid=$!

cleanup() {
  kill "${server_pid}" 2>/dev/null || true
}
trap cleanup EXIT

export SCHEMA_URL BASE_URL SCHEMATHESIS_STARTUP_TIMEOUT
ready=0
for _ in {1..120}; do
  if ! kill -0 "${server_pid}" 2>/dev/null; then
    echo "Schemathesis API server exited before becoming ready." >&2
    echo "Server log: ${SERVER_LOG}" >&2
    tail -n 80 "${SERVER_LOG}" >&2 || true
    break
  fi

  if "${PYTHON_BIN}" - <<'PY'
import os
import sys
import urllib.request

url = os.environ["SCHEMA_URL"]
timeout = float(os.environ["SCHEMATHESIS_STARTUP_TIMEOUT"])
try:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        sys.exit(0 if resp.status == 200 else 1)
except Exception:
    sys.exit(1)
PY
  then
    ready=1
    break
  fi
  sleep 0.25
done

if [ "${ready}" -ne 1 ]; then
  echo "Schema not available at ${SCHEMA_URL} after startup wait." >&2
  echo "Server log: ${SERVER_LOG}" >&2
  tail -n 80 "${SERVER_LOG}" >&2 || true
  exit 1
fi

max_examples="${SCHEMATHESIS_MAX_EXAMPLES:-}"
max_failures="${SCHEMATHESIS_MAX_FAILURES:-}"
phases="${SCHEMATHESIS_PHASES:-}"

case "${SCHEMATHESIS_MODE}" in
  smoke)
    max_examples="${max_examples:-5}"
    max_failures="${max_failures:-10}"
    phases="${phases:-coverage}"
    ;;
  deep)
    max_examples="${max_examples:-30}"
    max_failures="${max_failures:-100}"
    ;;
  *)
    echo "Invalid SCHEMATHESIS_MODE='${SCHEMATHESIS_MODE}' (expected: smoke or deep)" >&2
    exit 1
    ;;
esac

cmd=(
  schemathesis run "${SCHEMA_URL}"
  --url "${BASE_URL}"
  --checks "${SCHEMATHESIS_CHECKS:-not_a_server_error}"
  --include-path-regex '^/api/'
  --max-examples "${max_examples}"
  --max-failures "${max_failures}"
  --seed "${SCHEMATHESIS_SEED}"
  --generation-deterministic
)

if [ -n "${phases}" ]; then
  cmd+=(--phases "${phases}")
fi

if [ -n "${SCHEMATHESIS_AUTH_HEADER:-}" ]; then
  cmd+=(--header "${SCHEMATHESIS_AUTH_HEADER}")
elif [ -n "${SCHEMATHESIS_AUTH_TOKEN:-}" ]; then
  cmd+=(--header "Authorization: Bearer ${SCHEMATHESIS_AUTH_TOKEN}")
fi

if [ -n "${SCHEMATHESIS_TENANT_ID:-}" ]; then
  cmd+=(--header "x-tenant-id: ${SCHEMATHESIS_TENANT_ID}")
fi
if [ -n "${SCHEMATHESIS_DEV_USER_ID:-}" ]; then
  cmd+=(--header "X-Dev-User-Id: ${SCHEMATHESIS_DEV_USER_ID}")
fi
if [ -n "${SCHEMATHESIS_DEV_USER_EMAIL:-}" ]; then
  cmd+=(--header "X-Dev-User-Email: ${SCHEMATHESIS_DEV_USER_EMAIL}")
fi

if [ -n "${SCHEMATHESIS_EXTRA_ARGS:-}" ]; then
  # shellcheck disable=SC2206
  extra_args=( ${SCHEMATHESIS_EXTRA_ARGS} )
  cmd+=("${extra_args[@]}")
fi

echo "Running Schemathesis mode='${SCHEMATHESIS_MODE}' against ${SCHEMA_URL}" >&2
echo "Deterministic seed: ${SCHEMATHESIS_SEED}" >&2
echo "Server log: ${SERVER_LOG}" >&2
echo "To reproduce with identical data generation, reuse SCHEMATHESIS_SEED=${SCHEMATHESIS_SEED}" >&2

"${cmd[@]}"
