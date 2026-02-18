from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import main
from services.zpe import queue_targets as zpe_queue_targets


def _headers(user_id: str) -> dict[str, str]:
    return {
        "Authorization": f"Bearer {user_id}",
        "X-Tenant-Id": f"tenant-{user_id}",
    }


def _seed_target(
    store: zpe_queue_targets.QueueTargetStore,
    *,
    user_id: str,
    queue_name: str,
    server_id: str,
) -> str:
    target = store.add_target(
        user_id=user_id,
        queue_name=queue_name,
        server_id=server_id,
        name=f"{queue_name}-name",
    )
    return target.target_id


def test_queue_target_list_and_select(monkeypatch) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    store = zpe_queue_targets.QueueTargetStore(redis=fake)

    first_id = _seed_target(
        store,
        user_id="dev-user-1",
        queue_name="zpe",
        server_id="server-zpe",
    )
    second_id = _seed_target(
        store,
        user_id="dev-user-1",
        queue_name="zpe-2",
        server_id="server-zpe-2",
    )
    store.set_active_target("dev-user-1", first_id)

    client = TestClient(main.app)
    headers = _headers("dev-user-1")

    listing = client.get("/api/zpe/targets", headers=headers)
    assert listing.status_code == 200
    payload = listing.json()
    assert len(payload["targets"]) == 2
    assert payload["active_target_id"] == first_id

    select = client.put(f"/api/zpe/targets/{second_id}/active", headers=headers)
    assert select.status_code == 200
    assert select.json()["active_target_id"] == second_id

    listing2 = client.get("/api/zpe/targets", headers=headers)
    assert listing2.status_code == 200
    assert listing2.json()["active_target_id"] == second_id


def test_queue_target_is_user_scoped(monkeypatch) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    store = zpe_queue_targets.QueueTargetStore(redis=fake)

    owner_target = _seed_target(
        store,
        user_id="owner",
        queue_name="zpe-owner",
        server_id="server-owner",
    )
    store.set_active_target("owner", owner_target)

    client = TestClient(main.app)
    owner_headers = _headers("owner")
    other_headers = _headers("other")

    owner_listing = client.get("/api/zpe/targets", headers=owner_headers)
    assert owner_listing.status_code == 200
    assert len(owner_listing.json()["targets"]) == 1

    other_listing = client.get("/api/zpe/targets", headers=other_headers)
    assert other_listing.status_code == 200
    assert other_listing.json()["targets"] == []

    forbidden = client.put(f"/api/zpe/targets/{owner_target}/active", headers=other_headers)
    assert forbidden.status_code == 404
