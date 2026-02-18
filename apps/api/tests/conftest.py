from __future__ import annotations

from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import app.deps as deps  # noqa: E402
from services.authn.settings import AuthnSettings  # noqa: E402
from services.authn.types import UserIdentity  # noqa: E402
from services import structures as structure_service  # noqa: E402
from services.zpe.settings import ZPESettings  # noqa: E402


@pytest.fixture(autouse=True)
def _default_auth_mode_clerk(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        deps,
        "get_authn_settings",
        lambda: AuthnSettings(
            mode="clerk",
            clerk_issuer="https://issuer.example",
        ),
    )
    monkeypatch.setattr(
        deps,
        "verify_clerk_token",
        lambda token: UserIdentity(user_id=token, email=f"{token}@example.com"),
    )


@pytest.fixture(autouse=True)
def _default_structure_store_memory(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "memory")
    structure_service.reload_structure_store()


@pytest.fixture
def admin_auth_headers(monkeypatch: pytest.MonkeyPatch) -> dict[str, str]:
    token = "test-admin-token"
    monkeypatch.setattr(deps, "get_zpe_settings", lambda: ZPESettings(admin_token=token))
    return {"Authorization": f"Bearer {token}"}
