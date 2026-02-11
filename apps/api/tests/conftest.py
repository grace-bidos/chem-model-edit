from __future__ import annotations

from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import app.deps as deps  # noqa: E402
from services.authn.settings import AuthnSettings  # noqa: E402


@pytest.fixture(autouse=True)
def _default_auth_mode_dev_bypass(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        deps,
        "get_authn_settings",
        lambda: AuthnSettings(
            mode="dev-bypass",
            dev_bypass_user_id_header="X-Dev-User-Id",
            dev_bypass_email_header="X-Dev-User-Email",
        ),
    )
