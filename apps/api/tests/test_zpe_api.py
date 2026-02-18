from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

import main


LEGACY_REMOVED_ENDPOINTS = [
    ("post", "/api/zpe/jobs", {"calc_type": "qe.zpe.v1", "content": "", "mobile_indices": [0]}),
    ("get", "/api/zpe/jobs/job-1", None),
    ("get", "/api/zpe/jobs/job-1/result", None),
    ("get", "/api/zpe/jobs/job-1/files?kind=summary", None),
    ("post", "/api/zpe/compute/enroll-tokens", {"ttl_seconds": 60}),
    ("post", "/api/zpe/compute/servers", {"token": "x"}),
    ("delete", "/api/zpe/compute/servers/srv-1", None),
    ("post", "/api/zpe/compute/jobs/job-1/result", {"tenant_id": "t", "lease_id": "l", "result": {}, "summary_text": "", "freqs_csv": "", "meta": {}}),
    ("post", "/api/zpe/compute/jobs/job-1/failed", {"tenant_id": "t", "lease_id": "l", "error_code": "E", "error_message": "m", "traceback": "tb"}),
    ("post", "/api/zpe/compute/jobs/lease", None),
]


def test_zpe_parse_still_available() -> None:
    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/parse",
        json={
            "content": """
&control
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
 6.0 0.0 0.0
 0.0 6.0 0.0
 0.0 0.0 6.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0
 H 0.0 0.0 0.74
K_POINTS gamma
"""
        },
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["structure"]["atoms"][0]["symbol"] == "H"


@pytest.mark.parametrize(("method", "path", "body"), LEGACY_REMOVED_ENDPOINTS)
def test_legacy_compute_endpoints_removed(
    method: str,
    path: str,
    body: dict[str, object] | None,
) -> None:
    client = TestClient(main.app)
    requester = getattr(client, method)
    kwargs = {"json": body} if body is not None else {}
    response = requester(path, **kwargs)
    assert response.status_code == 404
