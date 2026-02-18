from __future__ import annotations

import pytest

from services.zpe import http_worker


def test_http_worker_is_retired_runtime_path() -> None:
    with pytest.raises(RuntimeError, match="retired"):
        http_worker.run_http_worker()
