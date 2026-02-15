from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from app import deps as app_deps
from services.zpe import enroll as zpe_enroll
from services.zpe import ops_flags as zpe_ops_flags
from services.zpe import result_store as zpe_store
from services.zpe import worker_auth as zpe_worker_auth
from services.zpe.settings import ZPESettings


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_ops_flags, "get_redis_connection", lambda: fake)
    return fake


def test_zpe_enroll_token_api(monkeypatch):
    _patch_redis(monkeypatch)
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(app_deps, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/enroll-tokens", json={"ttl_seconds": 60}
    )
    assert response.status_code == 401

    response = client.post(
        "/api/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60, "label": "lab"},
        headers={"Authorization": "Bearer secret"},
    )
    assert response.status_code == 200
    token = response.json()["token"]

    register = client.post(
        "/api/zpe/compute/servers",
        json={"token": token, "name": "server-1", "meta": {"gpu": 1}},
    )
    assert register.status_code == 200
    payload = register.json()
    assert payload["id"].startswith("compute-")
    assert payload["worker_token"]
    assert payload["token_expires_at"]
    assert payload["token_ttl_seconds"] > 0

    revoke = client.delete(
        f"/api/zpe/compute/servers/{payload['id']}",
        headers={"Authorization": "Bearer secret"},
    )
    assert revoke.status_code == 200
    assert revoke.json()["revoked_count"] == 1
