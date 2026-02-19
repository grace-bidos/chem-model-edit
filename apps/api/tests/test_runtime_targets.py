from __future__ import annotations

from fastapi.testclient import TestClient

import main
from services.runtime_nodes import RuntimeTarget


def _headers(user_id: str) -> dict[str, str]:
    return {
        "Authorization": f"Bearer {user_id}",
        "X-Tenant-Id": f"tenant-{user_id}",
    }


def _headers_for_tenant(user_id: str, tenant_id: str) -> dict[str, str]:
    return {
        "Authorization": f"Bearer {user_id}",
        "X-Tenant-Id": tenant_id,
    }


class _RuntimeNodeStore:
    def __init__(self) -> None:
        self._targets: dict[str, RuntimeTarget] = {}
        self._by_user: dict[tuple[str, str], list[str]] = {}
        self._active: dict[tuple[str, str], str] = {}

    def seed_target(self, *, tenant_id: str, user_id: str, queue_name: str, server_id: str) -> str:
        target = RuntimeTarget(
            target_id=f"qt-{len(self._targets) + 1}",
            tenant_id=tenant_id,
            user_id=user_id,
            queue_name=queue_name,
            server_id=server_id,
            registered_at="2026-02-18T00:00:00+00:00",
            name=f"{queue_name}-name",
        )
        self._targets[target.target_id] = target
        self._by_user.setdefault((tenant_id, user_id), []).append(target.target_id)
        return target.target_id

    def list_targets(self, tenant_id: str, user_id: str) -> list[RuntimeTarget]:
        ids = self._by_user.get((tenant_id, user_id), [])
        return [self._targets[target_id] for target_id in ids]

    def get_active_target(self, tenant_id: str, user_id: str) -> RuntimeTarget | None:
        target_id = self._active.get((tenant_id, user_id))
        if not target_id:
            return None
        return self._targets.get(target_id)

    def ensure_target_owner(self, tenant_id: str, user_id: str, target_id: str) -> RuntimeTarget:
        target = self._targets.get(target_id)
        if not target:
            raise KeyError("target not found")
        ids = self._by_user.get((tenant_id, user_id), [])
        if target_id not in ids:
            raise KeyError("target not found")
        return target

    def set_active_target(self, tenant_id: str, user_id: str, target_id: str) -> None:
        self._active[(tenant_id, user_id)] = target_id


def test_queue_target_list_and_select(monkeypatch) -> None:
    store = _RuntimeNodeStore()
    first_id = store.seed_target(
        tenant_id="tenant-dev-user-1",
        user_id="dev-user-1",
        queue_name="zpe",
        server_id="server-zpe",
    )
    second_id = store.seed_target(
        tenant_id="tenant-dev-user-1",
        user_id="dev-user-1",
        queue_name="zpe-2",
        server_id="server-zpe-2",
    )
    store.set_active_target("tenant-dev-user-1", "dev-user-1", first_id)

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: store)
    client = TestClient(main.app)
    headers = _headers("dev-user-1")

    listing = client.get("/api/runtime/targets", headers=headers)
    assert listing.status_code == 200
    payload = listing.json()
    assert len(payload["targets"]) == 2
    assert payload["active_target_id"] == first_id

    select = client.put(f"/api/runtime/targets/{second_id}/active", headers=headers)
    assert select.status_code == 200
    assert select.json()["active_target_id"] == second_id

    listing2 = client.get("/api/runtime/targets", headers=headers)
    assert listing2.status_code == 200
    assert listing2.json()["active_target_id"] == second_id


def test_queue_target_is_user_scoped(monkeypatch) -> None:
    store = _RuntimeNodeStore()
    owner_target = store.seed_target(
        tenant_id="tenant-owner",
        user_id="owner",
        queue_name="zpe-owner",
        server_id="server-owner",
    )
    store.set_active_target("tenant-owner", "owner", owner_target)

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: store)
    client = TestClient(main.app)
    owner_headers = _headers("owner")
    other_headers = _headers("other")

    owner_listing = client.get("/api/runtime/targets", headers=owner_headers)
    assert owner_listing.status_code == 200
    assert len(owner_listing.json()["targets"]) == 1

    other_listing = client.get("/api/runtime/targets", headers=other_headers)
    assert other_listing.status_code == 200
    assert other_listing.json()["targets"] == []

    forbidden = client.put(
        f"/api/runtime/targets/{owner_target}/active",
        headers=other_headers,
    )
    assert forbidden.status_code == 404


def test_queue_target_is_tenant_scoped_even_for_same_user_id(monkeypatch) -> None:
    store = _RuntimeNodeStore()
    tenant_a_target = store.seed_target(
        tenant_id="tenant-a",
        user_id="shared-user",
        queue_name="zpe-a",
        server_id="server-a",
    )
    store.set_active_target("tenant-a", "shared-user", tenant_a_target)

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: store)
    client = TestClient(main.app)
    tenant_a_headers = _headers_for_tenant("shared-user", "tenant-a")
    tenant_b_headers = _headers_for_tenant("shared-user", "tenant-b")

    tenant_a_listing = client.get("/api/runtime/targets", headers=tenant_a_headers)
    assert tenant_a_listing.status_code == 200
    assert len(tenant_a_listing.json()["targets"]) == 1

    tenant_b_listing = client.get("/api/runtime/targets", headers=tenant_b_headers)
    assert tenant_b_listing.status_code == 200
    assert tenant_b_listing.json()["targets"] == []

    forbidden = client.put(
        f"/api/runtime/targets/{tenant_a_target}/active",
        headers=tenant_b_headers,
    )
    assert forbidden.status_code == 404
