from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from services.zpe import enroll as zpe_enroll
from services.zpe import result_store as zpe_store
from services.zpe.settings import ZPESettings


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    return fake


def test_zpe_enroll_token_api(monkeypatch):
    _patch_redis(monkeypatch)
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(main, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    response = client.post("/calc/zpe/compute/enroll-tokens", json={"ttl_seconds": 60})
    assert response.status_code == 401

    response = client.post(
        "/calc/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60, "label": "lab"},
        headers={"Authorization": "Bearer secret"},
    )
    assert response.status_code == 200
    token = response.json()["token"]

    register = client.post(
        "/calc/zpe/compute/servers/register",
        json={"token": token, "name": "server-1", "meta": {"gpu": 1}},
    )
    assert register.status_code == 200
    assert register.json()["server_id"].startswith("compute-")
