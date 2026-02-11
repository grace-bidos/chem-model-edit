from __future__ import annotations

from fastapi.testclient import TestClient

import main


def test_legacy_auth_routes_are_removed() -> None:
    client = TestClient(main.app)

    assert client.post("/api/auth/register", json={}).status_code == 404
    assert client.post("/api/auth/login", json={}).status_code == 404
    assert client.post("/api/auth/logout").status_code == 404
    assert client.get("/api/auth/me").status_code == 404
