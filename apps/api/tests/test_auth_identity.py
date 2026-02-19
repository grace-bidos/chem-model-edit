from __future__ import annotations

from fastapi.testclient import TestClient

import app.deps as deps
import main
from services.authn.settings import AuthnSettings
from services.authn.types import UserIdentity

TENANT_HEADERS = {"X-Tenant-Id": "tenant-auth-tests"}


class _RuntimeNodeStore:
    def list_targets(self, _tenant_id: str, _user_id: str):  # type: ignore[no-untyped-def]
        return []

    def get_active_target(self, _tenant_id: str, _user_id: str):  # type: ignore[no-untyped-def]
        return None


def test_clerk_mode_accepts_verified_token(monkeypatch):
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
    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: _RuntimeNodeStore())

    client = TestClient(main.app)
    response = client.get(
        "/api/runtime/targets",
        headers={**TENANT_HEADERS, "Authorization": "Bearer valid-token"},
    )

    assert response.status_code == 200
    assert response.json()["targets"] == []


def test_clerk_mode_rejects_allowlist_denied(monkeypatch):
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
    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: _RuntimeNodeStore())

    client = TestClient(main.app)
    response = client.get(
        "/api/runtime/targets",
        headers={**TENANT_HEADERS, "Authorization": "Bearer denied-token"},
    )

    assert response.status_code == 403
