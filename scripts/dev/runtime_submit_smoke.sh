#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SECRETS_FILE="${DEV_SECRETS_FILE:-$ROOT_DIR/.tmp/dev-secrets.env}"
TENANT_ID="${TENANT_ID:-org_devtest}"

usage() {
  cat <<'USAGE'
Usage: runtime_submit_smoke.sh [--secrets-file <path>] [--tenant-id <tenant>]

Run a minimal submit smoke against FastAPI on Modal.
Required keys in secrets file:
- MODAL_API_BASE_URL
- MODAL_PROXY_KEY
- MODAL_PROXY_SECRET
- CLERK_TEST_JWT
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --secrets-file)
      shift
      [[ $# -gt 0 ]] || { echo "error: --secrets-file requires a value" >&2; exit 1; }
      SECRETS_FILE="$1"
      ;;
    --tenant-id)
      shift
      [[ $# -gt 0 ]] || { echo "error: --tenant-id requires a value" >&2; exit 1; }
      TENANT_ID="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
  shift
done

[[ -f "$SECRETS_FILE" ]] || { echo "error: secrets file not found: $SECRETS_FILE" >&2; exit 1; }

set -a
# shellcheck disable=SC1090
source "$SECRETS_FILE"
set +a

: "${MODAL_API_BASE_URL:?MODAL_API_BASE_URL is required in $SECRETS_FILE}"
: "${MODAL_PROXY_KEY:?MODAL_PROXY_KEY is required in $SECRETS_FILE}"
: "${MODAL_PROXY_SECRET:?MODAL_PROXY_SECRET is required in $SECRETS_FILE}"
: "${CLERK_TEST_JWT:?CLERK_TEST_JWT is required in $SECRETS_FILE}"

export TENANT_ID

cd "$ROOT_DIR/apps/api"
uv run python - <<'PY'
import base64
import json
import os
from urllib import error, request
import uuid

jwt = os.environ["CLERK_TEST_JWT"].strip()
parts = jwt.split(".")
if len(parts) < 2:
    raise SystemExit("invalid CLERK_TEST_JWT format")

payload_part = parts[1].replace("-", "+").replace("_", "/")
payload_part += "=" * ((4 - len(payload_part) % 4) % 4)
claims = json.loads(base64.b64decode(payload_part.encode("utf-8")))
sub = claims.get("sub")
if not isinstance(sub, str) or not sub:
    raise SystemExit("CLERK_TEST_JWT has no sub claim")

tenant = os.environ.get("TENANT_ID", "org_devtest")
job_id = f"job-dev-{uuid.uuid4().hex[:8]}"

payload = {
    "tenant_id": tenant,
    "workspace_id": tenant,
    "job_id": job_id,
    "idempotency_key": f"idem-{job_id}",
    "management_node_id": "mgmt-node-dev",
    "execution_profile": {"queue_name": "default"},
    "resource_shape": {"cpu": 1, "memory_mib": 512, "walltime_seconds": 600},
    "payload_ref": {"input_uri": f"inline://dev/{job_id}"},
    "requested_by": {"user_id": sub},
}

headers = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {jwt}",
    "x-tenant-id": tenant,
    "Modal-Key": os.environ["MODAL_PROXY_KEY"].strip(),
    "Modal-Secret": os.environ["MODAL_PROXY_SECRET"].strip(),
}

base = os.environ["MODAL_API_BASE_URL"].rstrip("/")
url = f"{base}/api/runtime/jobs:submit"
req = request.Request(url, data=json.dumps(payload).encode("utf-8"), headers=headers, method="POST")

try:
    with request.urlopen(req, timeout=20) as resp:
        body = json.loads(resp.read().decode("utf-8"))
        print(f"submit_status={resp.status}")
        print(f"job_id={body.get('job_id')}")
        print(f"submission_id={body.get('submission_id')}")
except error.HTTPError as exc:
    detail = exc.read().decode("utf-8", errors="replace")
    print(f"submit_status={exc.code}")
    print(detail)
    raise SystemExit(1)
PY
