from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_SUBMISSION_KEY = "zpe:ops:submission_enabled"
_DEQUEUE_KEY = "zpe:ops:dequeue_enabled"
_RESULT_READ_SOURCE_KEY = "zpe:ops:result_read_source"
ResultReadSource = Literal["redis", "projection"]


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


def _parse_result_read_source(
    raw: bytes | str | None, default: ResultReadSource
) -> ResultReadSource:
    if raw is None:
        return default
    if isinstance(raw, bytes):
        value = raw.decode("utf-8").strip().lower()
    else:
        value = str(raw).strip().lower()
    if value in {"redis", "projection"}:
        return cast(ResultReadSource, value)
    return default


@dataclass(frozen=True)
class OpsFlags:
    submission_enabled: bool
    dequeue_enabled: bool
    result_read_source: ResultReadSource


def _write_ops_value(
    *,
    redis: Redis,
    key: str,
    value: str,
    ttl_seconds: int,
) -> None:
    if ttl_seconds > 0:
        redis.setex(key, ttl_seconds, value)
        return
    redis.set(key, value)


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
    result_read_source = _parse_result_read_source(
        cast(bytes | str | None, redis.get(_RESULT_READ_SOURCE_KEY)),
        settings.cutover_result_read_source,
    )
    return OpsFlags(
        submission_enabled=submission_enabled,
        dequeue_enabled=dequeue_enabled,
        result_read_source=result_read_source,
    )


def set_ops_flags(
    *,
    submission_enabled: Optional[bool] = None,
    dequeue_enabled: Optional[bool] = None,
    result_read_source: Optional[ResultReadSource] = None,
    redis: Optional[Redis] = None,
) -> OpsFlags:
    settings = get_zpe_settings()
    ttl_seconds = settings.ops_flag_ttl_seconds
    redis = redis or get_redis_connection()
    if submission_enabled is not None:
        _write_ops_value(
            redis=redis,
            key=_SUBMISSION_KEY,
            value="1" if submission_enabled else "0",
            ttl_seconds=ttl_seconds,
        )
    if dequeue_enabled is not None:
        _write_ops_value(
            redis=redis,
            key=_DEQUEUE_KEY,
            value="1" if dequeue_enabled else "0",
            ttl_seconds=ttl_seconds,
        )
    if result_read_source is not None:
        _write_ops_value(
            redis=redis,
            key=_RESULT_READ_SOURCE_KEY,
            value=result_read_source,
            ttl_seconds=ttl_seconds,
        )
    return get_ops_flags(redis=redis)
