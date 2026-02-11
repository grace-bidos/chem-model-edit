from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
import app.routers.zpe as zpe_router
from services.zpe import backends as zpe_backends
from services.zpe import enroll as zpe_enroll
from services.zpe import job_meta as zpe_job_meta
from services.zpe import job_owner as zpe_job_owner
from services.zpe import ops_flags as zpe_ops_flags
from services.zpe import queue as zpe_queue
from services.zpe import queue_targets as zpe_queue_targets
from services.zpe import result_store as zpe_store
from services.zpe import worker_auth as zpe_worker_auth
from services.zpe.result_store import RedisResultStore
from services.zpe.settings import ZPESettings

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=1, ntyp=1
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
 H 0.0 0.0 0.0 1 1 1
""".strip()


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_owner, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_meta, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_ops_flags, "get_redis_connection", lambda: fake)
    return fake


def _setup_user_and_target(
    client: TestClient,
    monkeypatch,
    fake,
    *,
    user_id: str = "dev-user-1",
) -> tuple[dict[str, str], str]:
    headers = {
        "X-Dev-User-Id": user_id,
        "X-Dev-User-Email": f"{user_id}@example.com",
    }

    enroll = client.post(
        "/api/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60},
        headers=headers,
    )
    enroll_token = enroll.json()["token"]

    client.post(
        "/api/zpe/compute/servers",
        json={
            "token": enroll_token,
            "name": "server-1",
            "queue_name": "zpe",
        },
    )
    return headers, user_id


def test_zpe_mock_api_flow(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_router, "get_result_store", lambda: store)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    status = client.get(f"/api/zpe/jobs/{job_id}", headers=headers)
    assert status.status_code == 200
    assert status.json()["status"] == "finished"

    result = client.get(f"/api/zpe/jobs/{job_id}/result", headers=headers)
    assert result.status_code == 200
    payload = result.json()["result"]
    assert payload["zpe_ev"] >= 0.0
    assert payload["mobile_indices"] == [0]

    summary = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "summary"},
        headers=headers,
    )
    assert summary.status_code == 200
    assert "ZPE summary" in summary.text

    freqs = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "freqs"},
        headers=headers,
    )
    assert freqs.status_code == 200
    assert freqs.text.startswith("frequency_cm^-1")


def test_zpe_job_status_missing(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    monkeypatch.setattr(zpe_router, "get_result_store", lambda: store)
    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.get("/api/zpe/jobs/missing-job", headers=headers)
    assert response.status_code == 404


def test_zpe_job_result_not_finished(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    job_id = "job-not-finished"
    store.set_status(job_id, "queued")
    monkeypatch.setattr(zpe_router, "get_result_store", lambda: store)

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    owner_store = zpe_job_owner.JobOwnerStore(redis=fake)
    owner_store.set_owner(job_id, user_id)
    response = client.get(f"/api/zpe/jobs/{job_id}/result", headers=headers)
    assert response.status_code == 409


def test_zpe_submission_disabled(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_router, "get_result_store", lambda: store)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    zpe_ops_flags.set_ops_flags(submission_enabled=False, redis=fake)

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 503


def test_zpe_owner_enforcement_blocks_non_owner(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_router, "get_result_store", lambda: store)

    client = TestClient(main.app)
    owner_headers, _owner_id = _setup_user_and_target(
        client, monkeypatch, fake, user_id="user-owner"
    )
    other_headers = {
        "X-Dev-User-Id": "user-other",
        "X-Dev-User-Email": "user-other@example.com",
    }

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=owner_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    status_forbidden = client.get(f"/api/zpe/jobs/{job_id}", headers=other_headers)
    assert status_forbidden.status_code == 403

    result_forbidden = client.get(
        f"/api/zpe/jobs/{job_id}/result", headers=other_headers
    )
    assert result_forbidden.status_code == 403

    file_forbidden = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "summary"},
        headers=other_headers,
    )
    assert file_forbidden.status_code == 403
