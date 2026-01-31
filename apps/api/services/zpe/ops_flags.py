from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_SUBMISSION_KEY = "zpe:ops:submission_enabled"
_DEQUEUE_KEY = "zpe:ops:dequeue_enabled"


def _parse_bool(raw: bytes | str | None, default: bool) -> bool:
    if raw is None:
        return default
    if isinstance(raw, bytes):
        value = raw.decode("utf-8").strip().lower()
    else:
        value = str(raw).strip().lower()
    if value in {"1", "true", "yes", "on"}:
        return True
    if value in {"0", "false", "no", "off"}:
        return False
    return default


@dataclass(frozen=True)
class OpsFlags:
    submission_enabled: bool
    dequeue_enabled: bool


def get_ops_flags(*, redis: Optional[Redis] = None) -> OpsFlags:
    settings = get_zpe_settings()
    redis = redis or get_redis_connection()
    submission_enabled = _parse_bool(
        cast(bytes | str | None, redis.get(_SUBMISSION_KEY)),
        settings.submission_enabled,
    )
    dequeue_enabled = _parse_bool(
        cast(bytes | str | None, redis.get(_DEQUEUE_KEY)),
        settings.dequeue_enabled,
    )
    return OpsFlags(
        submission_enabled=submission_enabled,
        dequeue_enabled=dequeue_enabled,
    )


def set_ops_flags(
    *,
    submission_enabled: Optional[bool] = None,
    dequeue_enabled: Optional[bool] = None,
    redis: Optional[Redis] = None,
) -> OpsFlags:
    redis = redis or get_redis_connection()
    if submission_enabled is not None:
        redis.set(_SUBMISSION_KEY, "1" if submission_enabled else "0")
    if dequeue_enabled is not None:
        redis.set(_DEQUEUE_KEY, "1" if dequeue_enabled else "0")
    return get_ops_flags(redis=redis)
