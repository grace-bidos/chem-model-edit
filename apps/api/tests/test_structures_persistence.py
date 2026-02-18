from __future__ import annotations

from fastapi.testclient import TestClient

import main
from services import structures as structure_service

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0
 H 0.0 0.0 0.74
""".strip()


def test_memory_backend_persists_without_convex_credentials(
    monkeypatch,
) -> None:
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "memory")
    monkeypatch.delenv("CONVEX_URL", raising=False)
    monkeypatch.delenv("CONVEX_DEPLOY_KEY", raising=False)
    structure_service.reload_structure_store()

    client = TestClient(main.app)
    created = client.post("/api/structures", json={"content": QE_INPUT})
    assert created.status_code == 200
    structure_id = created.json()["structure_id"]

    fetched = client.get(f"/api/structures/{structure_id}")
    assert fetched.status_code == 200
    payload = fetched.json()
    assert payload["raw_input"] == QE_INPUT
    assert len(payload["structure"]["atoms"]) == 2


def test_convex_backend_rejects_missing_credentials(
    monkeypatch,
) -> None:
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "convex")
    monkeypatch.delenv("CONVEX_URL", raising=False)
    monkeypatch.delenv("CONVEX_DEPLOY_KEY", raising=False)
    try:
        try:
            structure_service.reload_structure_store()
        except RuntimeError as exc:
            assert "CONVEX_URL and CONVEX_DEPLOY_KEY are required" in str(exc)
        else:
            raise AssertionError("expected missing convex credentials to fail")
    finally:
        monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "memory")
        structure_service.reload_structure_store()
