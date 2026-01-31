#!/usr/bin/env bash
set -euo pipefail

HOST="127.0.0.1"
PORT="${SCHEMATHESIS_PORT:-8000}"
BASE_URL="http://${HOST}:${PORT}"
SCHEMA_URL="${BASE_URL}/openapi.json"

uvicorn main:app --host "${HOST}" --port "${PORT}" --log-level warning &
server_pid=$!

cleanup() {
  kill "${server_pid}" 2>/dev/null || true
}
trap cleanup EXIT

export SCHEMA_URL BASE_URL
ready=0
for _ in {1..40}; do
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

schemathesis run "${SCHEMA_URL}" --url "${BASE_URL}"
