from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from services.zpe import enroll as zpe_enroll
from services.zpe import job_owner as zpe_job_owner
from services.zpe import ops_flags as zpe_ops_flags
from services.zpe import queue as zpe_queue
from services.zpe import queue_targets as zpe_queue_targets
from services.zpe import result_store as zpe_store
from services.zpe import worker_auth as zpe_worker_auth


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_owner, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_ops_flags, "get_redis_connection", lambda: fake)
    return fake


def _auth_headers(client: TestClient, monkeypatch, fake) -> dict[str, str]:
    return {
        "X-Dev-User-Id": "dev-user-1",
        "X-Dev-User-Email": "dev-user@example.com",
        "X-Tenant-Id": "tenant-dev-user-1",
    }


def _register_target(
    client: TestClient,
    headers: dict[str, str],
    queue_name: str,
    activate_target: bool = False,
) -> str:
    enroll = client.post(
        "/api/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60},
        headers=headers,
    )
    enroll_token = enroll.json()["token"]
    register = client.post(
        "/api/zpe/compute/servers",
        json={
            "token": enroll_token,
            "name": f"server-{queue_name}",
            "queue_name": queue_name,
            "activate_target": activate_target,
        },
    )
    assert register.status_code == 200
    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.status_code == 200
    target = next(
        item for item in listing.json()["targets"] if item["queue_name"] == queue_name
    )
    return target["id"]


def test_queue_target_list_and_select(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    headers = _auth_headers(client, monkeypatch, fake)

    first_id = _register_target(client, headers, "zpe")

    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.status_code == 200
    payload = listing.json()
    assert len(payload["targets"]) == 1
    active_id = payload["active_target_id"]
    assert active_id is not None

    _register_target(client, headers, "zpe-2")

    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.status_code == 200
    payload = listing.json()
    assert len(payload["targets"]) == 2
    assert payload["active_target_id"] == first_id
    second = next(
        target for target in payload["targets"] if target["queue_name"] == "zpe-2"
    )
    second_id = second["id"]

    select = client.put(f"/api/zpe/targets/{second_id}/active", headers=headers)
    assert select.status_code == 200
    assert select.json()["active_target_id"] == second_id

    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.json()["active_target_id"] == second_id


def test_queue_target_visibility_is_scoped_by_user(monkeypatch):
    _patch_redis(monkeypatch)
    client = TestClient(main.app)
    owner_headers = {
        "X-Dev-User-Id": "user-owner",
        "X-Dev-User-Email": "user-owner@example.com",
        "X-Tenant-Id": "tenant-user-owner",
    }
    other_headers = {
        "X-Dev-User-Id": "user-other",
        "X-Dev-User-Email": "user-other@example.com",
        "X-Tenant-Id": "tenant-user-other",
    }

    owner_target_id = _register_target(client, owner_headers, "zpe-owner")

    owner_listing = client.get("/api/zpe/targets", headers=owner_headers)
    assert owner_listing.status_code == 200
    assert len(owner_listing.json()["targets"]) == 1

    other_listing = client.get("/api/zpe/targets", headers=other_headers)
    assert other_listing.status_code == 200
    assert other_listing.json()["targets"] == []

    forbidden_select = client.put(
        f"/api/zpe/targets/{owner_target_id}/active", headers=other_headers
    )
    assert forbidden_select.status_code == 404


def test_register_target_activate_flag_can_switch_active(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    headers = _auth_headers(client, monkeypatch, fake)

    first_id = _register_target(client, headers, "zpe")
    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.json()["active_target_id"] == first_id

    second_id = _register_target(client, headers, "zpe-2", activate_target=True)
    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.json()["active_target_id"] == second_id
