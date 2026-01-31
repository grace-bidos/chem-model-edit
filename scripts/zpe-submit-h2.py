#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path


def main() -> None:
    base_url = os.getenv("ZPE_API_URL", "http://localhost:8000").rstrip("/")
    if not base_url.endswith("/api"):
        base_url = f"{base_url}/api"
    content = Path("samples/qe-in/h2_zpe.in").read_text()
    payload = {
        "content": content,
        "mobileIndices": [0, 1],
        "calcMode": "new",
        "useEnviron": False,
    }
    req = urllib.request.Request(
        f"{base_url}/zpe/jobs",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req) as resp:
        print(resp.read().decode("utf-8"))


if __name__ == "__main__":
    main()
