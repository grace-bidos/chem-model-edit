from __future__ import annotations

from datetime import datetime, timezone
from io import BytesIO
import json
from urllib import error as urlerror

from fastapi.testclient import TestClient

import main
from services import structures as structure_service
from app.schemas.common import QeParameters

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


def test_inmemory_store_tenant_workspace_scope_filter() -> None:
    store = structure_service.InMemoryStructureStore()
    structure_id = store.create(
        atoms=structure_service.parse_qe_atoms(QE_INPUT)[0],
        source="qe",
        cif="data_test",
        params=None,
        raw_input=QE_INPUT,
        tenant_id="t1",
        workspace_id="w1",
    )
    assert store.get(structure_id, tenant_id="t1", workspace_id="w1") is not None
    assert store.get(structure_id, tenant_id="t2") is None
    assert store.get(structure_id, workspace_id="w2") is None


def test_convex_request_payload_variants(monkeypatch) -> None:
    store = structure_service.ConvexStructureStore(
        convex_url="https://convex.example",
        convex_deploy_key="dev:key",
        timeout_seconds=3,
    )

    class _Resp:
        def __init__(self, payload: str) -> None:
            self._payload = payload.encode("utf-8")

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self) -> bytes:
            return self._payload

    responses = iter(
        [
            _Resp("{}"),
            _Resp(json.dumps({"value": [1, 2]})),
            _Resp(json.dumps({"value": {"ok": True}})),
        ]
    )
    monkeypatch.setattr(structure_service.urlrequest, "urlopen", lambda *_a, **_k: next(responses))

    assert store._request("query", "x:y", {}) == {}
    assert store._request("query", "x:y", {}) == {"value": [1, 2]}
    assert store._request("query", "x:y", {}) == {"ok": True}


def test_convex_request_errors_and_invalid_payload(monkeypatch) -> None:
    store = structure_service.ConvexStructureStore(
        convex_url="https://convex.example",
        convex_deploy_key="dev:key",
    )

    def raise_http(*_a, **_k):
        raise urlerror.HTTPError(
            url="https://convex.example/api/query",
            code=500,
            msg="boom",
            hdrs=None,
            fp=BytesIO(b"bad"),
        )

    monkeypatch.setattr(structure_service.urlrequest, "urlopen", raise_http)
    try:
        store._request("query", "x:y", {})
    except RuntimeError as exc:
        assert "HTTP 500" in str(exc)
    else:
        raise AssertionError("expected http error")

    monkeypatch.setattr(
        structure_service.urlrequest,
        "urlopen",
        lambda *_a, **_k: (_ for _ in ()).throw(urlerror.URLError("down")),
    )
    try:
        store._request("query", "x:y", {})
    except RuntimeError as exc:
        assert "convex query failed" in str(exc)
    else:
        raise AssertionError("expected url error")

    class _Resp:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self) -> bytes:
            return b"[]"

    monkeypatch.setattr(structure_service.urlrequest, "urlopen", lambda *_a, **_k: _Resp())
    try:
        store._request("query", "x:y", {})
    except RuntimeError as exc:
        assert "invalid payload" in str(exc)
    else:
        raise AssertionError("expected invalid payload")


def test_convex_create_and_get_roundtrip_logic(monkeypatch) -> None:
    store = structure_service.ConvexStructureStore(
        convex_url="https://convex.example",
        convex_deploy_key="dev:key",
        chunk_size=10,
    )
    atoms = structure_service.parse_qe_atoms(QE_INPUT)[0]

    mutations: list[tuple[str, dict]] = []
    monkeypatch.setattr(store, "_mutation", lambda path, args: mutations.append((path, args)) or {})
    monkeypatch.setattr(
        structure_service,
        "_atoms_from_cif",
        lambda _cif: atoms,
    )
    structure_id = store.create(
        atoms=atoms,
        source="qe",
        cif="data_test",
        params=QeParameters(),
        raw_input="1234567890ABCDEFG",
        tenant_id="tenant-a",
        workspace_id="ws-a",
    )
    assert structure_id
    assert mutations[0][0] == "structures:putStructureRecord"
    assert mutations[1][0] == "structures:putStructureRawChunk"

    queries = {
        "structures:getStructureRecordById": {
            "_id": "doc1",
            "tenant_id": "tenant-a",
            "workspace_id": "ws-a",
            "source": "qe",
            "cif": "data_test",
            "params": {"control": {}, "system": {}, "electrons": {}, "ions": {}, "cell": {}},
            "created_at": "2026-02-18T00:00:00",
        },
        "structures:getStructureRawChunks": {
            "value": [
                {"chunk_index": 1, "chunk_text": "DEF"},
                {"chunk_index": 0, "chunk_text": "ABC"},
            ]
        },
    }
    monkeypatch.setattr(store, "_query", lambda path, _args: queries[path])
    entry = store.get(structure_id, tenant_id="tenant-a", workspace_id="ws-a")
    assert entry is not None
    assert entry.raw_input == "ABCDEF"
    assert entry.params is not None
    assert entry.created_at.tzinfo == timezone.utc

    queries["structures:getStructureRawChunks"] = [
        {"chunk_index": 0, "chunk_text": "X"},
        {"chunk_index": 1, "chunk_text": "Y"},
    ]
    assert store.get(structure_id) is not None
    queries["structures:getStructureRecordById"] = {"_id": None}
    assert store.get(structure_id) is None


def test_convex_get_scope_mismatch_returns_none(monkeypatch) -> None:
    store = structure_service.ConvexStructureStore(
        convex_url="https://convex.example",
        convex_deploy_key="dev:key",
    )
    monkeypatch.setattr(
        store,
        "_query",
        lambda _p, _a: {
            "_id": "doc1",
            "tenant_id": "tenant-a",
            "workspace_id": "ws-a",
            "source": "qe",
            "cif": "data_test",
            "created_at": datetime.now(timezone.utc).isoformat(),
        },
    )
    assert store.get("sid", tenant_id="tenant-b") is None
    assert store.get("sid", workspace_id="ws-b") is None


def test_build_store_configuration_paths(monkeypatch) -> None:
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "invalid")
    try:
        structure_service._build_store()
    except ValueError as exc:
        assert "STRUCTURE_STORE_BACKEND" in str(exc)
    else:
        raise AssertionError("expected invalid backend error")

    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "convex")
    monkeypatch.setenv("CONVEX_URL", "https://convex.example")
    monkeypatch.setenv("CONVEX_DEPLOY_KEY", "dev:key")
    monkeypatch.setenv("STRUCTURE_RAW_CHUNK_SIZE", "100")
    monkeypatch.setenv("STRUCTURE_STORE_TIMEOUT_SECONDS", "7")
    built = structure_service._build_store()
    assert isinstance(built, structure_service.ConvexStructureStore)


def test_wrapper_accessors_and_register(monkeypatch) -> None:
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "memory")
    structure_service.reload_structure_store()
    sid, _structure, _source, _params = structure_service.create_structure_from_qe(QE_INPUT)

    assert structure_service.get_structure_cif(sid)
    assert structure_service.get_structure_entry(sid)
    assert structure_service.get_structure(sid).atoms
    assert structure_service.get_structure_params(sid) is not None
    assert structure_service.get_structure_raw_input(sid) is not None

    atoms = structure_service.parse_qe_atoms(QE_INPUT)[0]
    sid2 = structure_service.register_structure_atoms(atoms, "qe")
    assert structure_service.get_structure_entry(sid2).raw_input is None

    try:
        structure_service.get_structure_entry("missing-id")
    except KeyError:
        pass
    else:
        raise AssertionError("expected missing-id error")
