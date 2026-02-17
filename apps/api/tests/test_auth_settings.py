from __future__ import annotations

from types import SimpleNamespace

import pytest

from services.auth import settings


@pytest.fixture(autouse=True)
def _clear_auth_settings_cache() -> None:
    settings.get_auth_settings.cache_clear()
    yield
    settings.get_auth_settings.cache_clear()


def test_build_auth_settings_uses_resolved_env_file(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
    env_path = "/tmp/auth-settings.env"

    class SpyAuthSettings(settings.AuthSettings):
        def __init__(self, **kwargs):
            captured.update(kwargs)
            super().__init__(**kwargs)

    monkeypatch.setattr(settings, "resolve_env_file", lambda: env_path)
    monkeypatch.setattr(settings, "AuthSettings", SpyAuthSettings)

    built = settings._build_auth_settings()

    assert isinstance(built, settings.AuthSettings)
    assert captured["_env_file"] == env_path


def test_get_auth_settings_falls_back_to_zpe_redis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        settings,
        "_build_auth_settings",
        lambda: settings.AuthSettings(redis_url=None),
    )
    monkeypatch.setattr(
        settings,
        "get_zpe_settings",
        lambda: SimpleNamespace(redis_url="redis://zpe.example:6379/9"),
    )

    resolved = settings.get_auth_settings()

    assert resolved.redis_url == "redis://zpe.example:6379/9"


def test_get_auth_settings_keeps_auth_redis_and_skips_zpe_lookup(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        settings,
        "_build_auth_settings",
        lambda: settings.AuthSettings(redis_url="redis://auth.example:6379/4"),
    )

    def _should_not_be_called():
        raise AssertionError("zpe settings lookup should not run when auth redis exists")

    monkeypatch.setattr(settings, "get_zpe_settings", _should_not_be_called)

    resolved = settings.get_auth_settings()

    assert resolved.redis_url == "redis://auth.example:6379/4"
