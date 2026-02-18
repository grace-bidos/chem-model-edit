from __future__ import annotations

import fakeredis
from fastapi.testclient import TestClient

import app.deps as deps
import main
from services.authn.settings import AuthnSettings
from services.authn.types import UserIdentity
from services.zpe import queue_targets as zpe_queue_targets

TENANT_HEADERS = {"X-Tenant-Id": "tenant-auth-tests"}


def _patch_target_store_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    return fake


def test_clerk_mode_accepts_verified_token(monkeypatch):
    _patch_target_store_redis(monkeypatch)
    monkeypatch.setattr(
        deps,
        "get_authn_settings",
        lambda: AuthnSettings(
            mode="clerk",
            clerk_issuer="https://issuer.example",
            allowed_emails="user@example.com",
        ),
    )
    monkeypatch.setattr(
        deps,
        "verify_clerk_token",
        lambda token: UserIdentity(user_id="user-clerk-1", email="user@example.com"),
    )

    client = TestClient(main.app)
    response = client.get(
        "/api/zpe/targets",
        headers={**TENANT_HEADERS, "Authorization": "Bearer valid-token"},
    )

    assert response.status_code == 200
    assert response.json()["targets"] == []


def test_clerk_mode_rejects_allowlist_denied(monkeypatch):
    _patch_target_store_redis(monkeypatch)
    monkeypatch.setattr(
        deps,
        "get_authn_settings",
        lambda: AuthnSettings(
            mode="clerk",
            clerk_issuer="https://issuer.example",
            allowed_emails="allowed@example.com",
        ),
    )

    def _deny(_token: str) -> UserIdentity:
        raise PermissionError("allowlist denied")

    monkeypatch.setattr(deps, "verify_clerk_token", _deny)

    client = TestClient(main.app)
    response = client.get(
        "/api/zpe/targets",
        headers={**TENANT_HEADERS, "Authorization": "Bearer denied-token"},
    )

    assert response.status_code == 403
