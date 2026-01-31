from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def new_request_id() -> str:
    return uuid4().hex


def log_event(
    logger: logging.Logger,
    *,
    event: str,
    level: int = logging.INFO,
    **fields: Any,
) -> None:
    payload = {"ts": _now_iso(), "event": event, **fields}
    logger.log(
        level,
        json.dumps(payload, ensure_ascii=True, separators=(",", ":")),
    )
