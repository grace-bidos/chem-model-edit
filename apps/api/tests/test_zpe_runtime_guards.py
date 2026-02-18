from __future__ import annotations

from types import SimpleNamespace

import pytest

from services.zpe import backends
from services.zpe import worker


def test_enqueue_zpe_job_rejects_unsupported_compute_mode(monkeypatch) -> None:
    monkeypatch.setattr(
        backends,
        "get_zpe_settings",
        lambda: SimpleNamespace(compute_mode="remote-queue"),
    )

    with pytest.raises(ValueError, match="compute_mode must be 'remote-http' or 'mock'"):
        backends.enqueue_zpe_job({"content": "", "mobile_indices": [], "calc_type": "zpe"})


def test_get_current_job_returns_none() -> None:
    assert worker.get_current_job() is None
