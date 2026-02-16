#!/usr/bin/env bash
set -euo pipefail

HOST="127.0.0.1"
if [ -n "${SCHEMATHESIS_PORT:-}" ]; then
  PORT="${SCHEMATHESIS_PORT}"
else
  PORT="$(
    python - <<'PY'
import socket

with socket.socket() as sock:
    sock.bind(("127.0.0.1", 0))
    print(sock.getsockname()[1])
PY
  )"
fi
BASE_URL="http://${HOST}:${PORT}"
SCHEMA_URL="${BASE_URL}/api/openapi.json"
MODE="${SCHEMATHESIS_MODE:-smoke}"
CHECKS="${SCHEMATHESIS_CHECKS:-not_a_server_error}"
MAX_EXAMPLES="${SCHEMATHESIS_MAX_EXAMPLES:-30}"
MAX_FAILURES="${SCHEMATHESIS_MAX_FAILURES:-5}"
SEED="${SCHEMATHESIS_SEED:-137}"
TENANT_ID="${SCHEMATHESIS_TENANT_ID:-tenant-dev}"
AUTH_TOKEN="${SCHEMATHESIS_AUTH_TOKEN:-}"

uvicorn main:app --host "${HOST}" --port "${PORT}" --log-level warning &
server_pid=$!

cleanup() {
  kill "${server_pid}" 2>/dev/null || true
}
trap cleanup EXIT

export SCHEMA_URL BASE_URL
ready=0
for _ in {1..40}; do
  if ! kill -0 "${server_pid}" 2>/dev/null; then
    echo "Schemathesis API server exited before becoming ready" >&2
    break
  fi

  if python - <<'PY'
import os
import sys
import urllib.request

url = os.environ["SCHEMA_URL"]
try:
    with urllib.request.urlopen(url) as resp:
        sys.exit(0 if resp.status == 200 else 1)
except Exception:
    sys.exit(1)
PY
  then
    ready=1
    break
  fi
  sleep 0.5
done

if [ "$ready" -ne 1 ]; then
  echo "Schema not available at ${SCHEMA_URL}" >&2
  exit 1
fi

SMOKE_REGEX='^/api/(health|structures/parse|structures|structures/[^/]+|structures/[^/]+/view|structures/export|transforms/delta-transplant|supercells/builds)$'
BROAD_REGEX='^/api/(health|structures.*|transforms.*|supercells.*|zpe.*)$'

if [ "${MODE}" = "broad" ]; then
  INCLUDE_REGEX="${BROAD_REGEX}"
else
  INCLUDE_REGEX="${SMOKE_REGEX}"
fi

cmd=(
  schemathesis run "${SCHEMA_URL}"
  --url "${BASE_URL}"
  --checks "${CHECKS}"
  --include-path-regex "${INCLUDE_REGEX}"
  --seed "${SEED}"
  --max-examples "${MAX_EXAMPLES}"
  --max-failures "${MAX_FAILURES}"
  --header "x-tenant-id: ${TENANT_ID}"
)

if [ -n "${AUTH_TOKEN}" ]; then
  cmd+=(--header "authorization: Bearer ${AUTH_TOKEN}")
fi

if [ -n "${SCHEMATHESIS_EXTRA_ARGS:-}" ]; then
  # shellcheck disable=SC2206
  extra_args=( ${SCHEMATHESIS_EXTRA_ARGS} )
  cmd+=("${extra_args[@]}")
fi

echo "Schemathesis mode=${MODE} include=${INCLUDE_REGEX} seed=${SEED}"
"${cmd[@]}"
