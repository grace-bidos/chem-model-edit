from __future__ import annotations

import json

import fakeredis
from fastapi.testclient import TestClient

import main
import app.routers.zpe as zpe_router
from services.zpe import lease as zpe_lease
from services.zpe import job_meta as zpe_job_meta
from services.zpe import ops_flags as zpe_ops_flags
from services.zpe import result_store as zpe_store
from services.zpe import worker_auth as zpe_worker_auth
from services.zpe.settings import ZPESettings


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_lease, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_meta, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_ops_flags, "get_redis_connection", lambda: fake)
    return fake


def test_lease_endpoint_returns_job(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(lease_ttl_seconds=60, result_ttl_seconds=60)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_lease, "get_zpe_settings", lambda: settings)

    token_store = zpe_worker_auth.WorkerTokenStore(redis=fake)
    worker_token = token_store.create_token("compute-1").token

    payload = {"content": "test", "mobile_indices": [0]}
    job_id = "http-test"
    fake.setex(f"zpe:payload:{job_id}", 60, json.dumps(payload))
    fake.lpush("zpe:queue", job_id)
    zpe_store.RedisResultStore(redis=fake).set_status(job_id, "queued")

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/jobs/lease",
        headers={"Authorization": f"Bearer {worker_token}"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["job_id"] == job_id
    assert data["payload"]["content"] == payload["content"]
    assert data["payload"]["mobile_indices"] == payload["mobile_indices"]
    assert data["lease_id"]
    assert data["lease_ttl_seconds"] == 60


def test_lease_endpoint_no_job(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(lease_ttl_seconds=60, result_ttl_seconds=60)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_lease, "get_zpe_settings", lambda: settings)

    token_store = zpe_worker_auth.WorkerTokenStore(redis=fake)
    worker_token = token_store.create_token("compute-1").token

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/jobs/lease",
        headers={"Authorization": f"Bearer {worker_token}"},
    )
    assert response.status_code == 204


def test_lease_endpoint_paused(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(lease_ttl_seconds=60, result_ttl_seconds=60)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_lease, "get_zpe_settings", lambda: settings)

    token_store = zpe_worker_auth.WorkerTokenStore(redis=fake)
    worker_token = token_store.create_token("compute-1").token
    zpe_ops_flags.set_ops_flags(dequeue_enabled=False, redis=fake)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/jobs/lease",
        headers={"Authorization": f"Bearer {worker_token}"},
    )
    assert response.status_code == 204
