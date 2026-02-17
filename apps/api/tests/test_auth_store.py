from __future__ import annotations

import json

import fakeredis
import pytest

from services.auth.settings import AuthSettings
from services.auth.store import AuthStore
import services.auth.store as auth_store


def _auth_settings(*, pepper: str | None = "pepper-1", ttl: int = 600) -> AuthSettings:
    return AuthSettings(
        redis_url="redis://unused:6379/0",
        password_iterations=2_000,
        password_pepper=pepper,
        session_ttl_seconds=ttl,
    )


def test_create_user_success_persists_index_and_user_payload(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(auth_store, "get_auth_settings", lambda: _auth_settings())
    store = AuthStore(redis=fake)

    user = store.create_user("  Test.User@Example.com  ", "strong-pass-123")

    assert user.user_id.startswith("user-")
    assert user.email == "test.user@example.com"
    assert user.password_hash.startswith("pbkdf2_sha256$")
    assert fake.get("auth:user_email:test.user@example.com") == user.user_id.encode(
        "utf-8"
    )

    raw_user = fake.get(f"auth:user:{user.user_id}")
    assert raw_user is not None
    payload = json.loads(raw_user.decode("utf-8"))
    assert payload["email"] == "test.user@example.com"
    assert payload["user_id"] == user.user_id


def test_create_user_rejects_duplicate_email_case_insensitively(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(auth_store, "get_auth_settings", lambda: _auth_settings())
    store = AuthStore(redis=fake)

    store.create_user("Person@Example.com", "first-pass-123")

    with pytest.raises(ValueError, match="email already registered"):
        store.create_user(" person@example.com ", "second-pass-456")


def test_authenticate_success_and_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(auth_store, "get_auth_settings", lambda: _auth_settings())
    store = AuthStore(redis=fake)
    created = store.create_user("auth@example.com", "valid-password")

    authenticated = store.authenticate(" AUTH@example.com ", "valid-password")
    bad_password = store.authenticate("auth@example.com", "invalid-password")
    missing_user = store.authenticate("missing@example.com", "any-password")

    assert authenticated is not None
    assert authenticated.user_id == created.user_id
    assert bad_password is None
    assert missing_user is None


def test_session_create_get_delete_lifecycle_with_fakeredis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(auth_store, "get_auth_settings", lambda: _auth_settings(ttl=90))
    store = AuthStore(redis=fake)

    session = store.create_session("user-abc123")
    fetched = store.get_user_by_session(session.token, refresh=False)
    store.delete_session(session.token)
    deleted = store.get_user_by_session(session.token, refresh=False)

    assert session.token
    assert session.user_id == "user-abc123"
    assert fetched is not None
    assert fetched.user_id == "user-abc123"
    assert fetched.token == session.token
    assert deleted is None
