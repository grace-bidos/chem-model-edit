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


def test_lease_endpoint_legacy_worker_disabled(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(lease_ttl_seconds=60, result_ttl_seconds=60)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_lease, "get_zpe_settings", lambda: settings)

    token_store = zpe_worker_auth.WorkerTokenStore(redis=fake)
    worker_token = token_store.create_token("compute-1").token
    zpe_ops_flags.set_ops_flags(legacy_worker_endpoints_enabled=False, redis=fake)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/jobs/lease",
        headers={"Authorization": f"Bearer {worker_token}"},
    )
    assert response.status_code == 503


def test_lease_endpoint_skips_out_of_scope_job_for_scoped_worker(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(lease_ttl_seconds=60, result_ttl_seconds=60)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_lease, "get_zpe_settings", lambda: settings)

    token_store = zpe_worker_auth.WorkerTokenStore(redis=fake)
    worker_token = token_store.create_token(
        "compute-1",
        tenant_id="tenant-1",
        workspace_id="workspace-1",
    ).token

    out_of_scope_job = "job-out-of-scope"
    out_of_scope_payload = {"content": "out", "mobile_indices": [0]}
    fake.setex(f"zpe:payload:{out_of_scope_job}", 60, json.dumps(out_of_scope_payload))
    zpe_store.RedisResultStore(redis=fake).set_status(out_of_scope_job, "queued")
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        out_of_scope_job,
        {
            "tenant_id": "tenant-2",
            "workspace_id": "workspace-2",
            "request_id": "req-out",
            "user_id": "user-out",
        },
    )

    in_scope_job = "job-in-scope"
    in_scope_payload = {"content": "in", "mobile_indices": [1]}
    fake.setex(f"zpe:payload:{in_scope_job}", 60, json.dumps(in_scope_payload))
    zpe_store.RedisResultStore(redis=fake).set_status(in_scope_job, "queued")
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        in_scope_job,
        {
            "tenant_id": "tenant-1",
            "workspace_id": "workspace-1",
            "request_id": "req-in",
            "user_id": "user-in",
        },
    )

    # RPOP consumes the rightmost item first, so enqueue out-of-scope first.
    fake.lpush("zpe:queue", out_of_scope_job)
    fake.lpush("zpe:queue", in_scope_job)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/compute/jobs/lease",
        headers={"Authorization": f"Bearer {worker_token}"},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["job_id"] == in_scope_job
    assert data["payload"]["content"] == in_scope_payload["content"]

    # Out-of-scope job is requeued and lease state is released.
    assert fake.hgetall(f"zpe:lease:{out_of_scope_job}") == {}
    assert fake.zscore("zpe:lease:index", out_of_scope_job) is None
    assert fake.llen("zpe:queue") == 1
