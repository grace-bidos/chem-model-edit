#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path
from uuid import uuid4


def main() -> None:
    base_url = os.getenv("ZPE_API_URL", "http://localhost:8000").rstrip("/")
    if not base_url.endswith("/api"):
        base_url = f"{base_url}/api"
    _ = Path("samples/qe-in/h2_zpe.in").read_text()
    tenant_id = os.getenv("ZPE_TENANT_ID", "tenant-dev")
    user_id = os.getenv("ZPE_USER_ID", "user-dev")
    auth_token = os.getenv("ZPE_AUTH_TOKEN", "")
    management_node_id = os.getenv("ZPE_MANAGEMENT_NODE_ID", "mgmt-dev")
    queue_name = os.getenv("ZPE_QUEUE_NAME", "default")
    job_id = f"job-{uuid4()}"
    payload = {
        "tenant_id": tenant_id,
        "workspace_id": tenant_id,
        "job_id": job_id,
        "idempotency_key": f"submit-{job_id}",
        "management_node_id": management_node_id,
        "execution_profile": {"queue_name": queue_name},
        "resource_shape": {"cpu": 1, "memory_mib": 2048, "walltime_seconds": 3600},
        "payload_ref": {"input_uri": f"inline://samples/h2/{job_id}"},
        "requested_by": {"user_id": user_id},
    }
    req = urllib.request.Request(
        f"{base_url}/runtime/jobs:submit",
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "X-Tenant-Id": tenant_id,
            **({"Authorization": f"Bearer {auth_token}"} if auth_token else {}),
        },
    )
    with urllib.request.urlopen(req) as resp:
        print(resp.read().decode("utf-8"))


if __name__ == "__main__":
    main()
