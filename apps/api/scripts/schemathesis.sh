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
SMOKE_CHECKS="not_a_server_error,content_type_conformance,response_headers_conformance"
BROAD_CHECKS="not_a_server_error,status_code_conformance,content_type_conformance,response_headers_conformance,response_schema_conformance,negative_data_rejection"
if [ "${MODE}" = "broad" ]; then
  DEFAULT_CHECKS="${BROAD_CHECKS}"
  DEFAULT_MAX_EXAMPLES=30
else
  DEFAULT_CHECKS="${SMOKE_CHECKS}"
  DEFAULT_MAX_EXAMPLES=8
fi
CHECKS="${SCHEMATHESIS_CHECKS:-${DEFAULT_CHECKS}}"
MAX_EXAMPLES="${SCHEMATHESIS_MAX_EXAMPLES:-${DEFAULT_MAX_EXAMPLES}}"
MAX_FAILURES="${SCHEMATHESIS_MAX_FAILURES:-5}"
SEED="${SCHEMATHESIS_SEED:-137}"
TENANT_ID="${SCHEMATHESIS_TENANT_ID:-tenant-dev}"
AUTH_TOKEN="${SCHEMATHESIS_AUTH_TOKEN:-}"
DETERMINISTIC="${SCHEMATHESIS_DETERMINISTIC:-1}"
GEN_DB="${SCHEMATHESIS_GENERATION_DB:-.schemathesis/examples.db}"
REPORT_DIR="${SCHEMATHESIS_REPORT_DIR:-schemathesis-report/${MODE}}"

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

if [ "${DETERMINISTIC}" = "1" ]; then
  cmd+=(--generation-deterministic)
  effective_gen_db="disabled"
else
  cmd+=(--generation-database "${GEN_DB}")
  effective_gen_db="${GEN_DB}"
fi

if [ -n "${SCHEMATHESIS_REPORTS:-}" ]; then
  cmd+=(--report "${SCHEMATHESIS_REPORTS}")
  cmd+=(--report-dir "${REPORT_DIR}")
fi

if [ -n "${AUTH_TOKEN}" ]; then
  cmd+=(--header "authorization: Bearer ${AUTH_TOKEN}")
fi

if [ -n "${SCHEMATHESIS_EXTRA_ARGS:-}" ]; then
  # shellcheck disable=SC2206
  extra_args=( ${SCHEMATHESIS_EXTRA_ARGS} )
  cmd+=("${extra_args[@]}")
fi

schemathesis_version="$(schemathesis --version | tr '\n' ' ' | sed 's/[[:space:]]\+/ /g' | sed 's/[[:space:]]$//')"
echo "Schemathesis mode=${MODE} include=${INCLUDE_REGEX} checks=${CHECKS} seed=${SEED} max_examples=${MAX_EXAMPLES} max_failures=${MAX_FAILURES} deterministic=${DETERMINISTIC} generation_db=${effective_gen_db} reports=${SCHEMATHESIS_REPORTS:-none} report_dir=${REPORT_DIR} tenant=${TENANT_ID} auth_header=$([ -n "${AUTH_TOKEN}" ] && echo present || echo absent) version=${schemathesis_version}"
echo "Schemathesis command: ${cmd[*]}"
"${cmd[@]}"
