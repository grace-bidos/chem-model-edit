from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_SUBMISSION_KEY = "zpe:ops:submission_enabled"
_DEQUEUE_KEY = "zpe:ops:dequeue_enabled"
_SUBMISSION_ROUTE_KEY = "zpe:ops:submission_route"
_RESULT_READ_SOURCE_KEY = "zpe:ops:result_read_source"
_LEGACY_WORKER_ENDPOINTS_KEY = "zpe:ops:legacy_worker_endpoints_enabled"

SubmissionRoute = Literal["redis-worker", "next-gen"]
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


def _parse_submission_route(
    raw: bytes | str | None, default: SubmissionRoute
) -> SubmissionRoute:
    if raw is None:
        return default
    if isinstance(raw, bytes):
        value = raw.decode("utf-8").strip().lower()
    else:
        value = str(raw).strip().lower()
    if value in {"redis-worker", "next-gen"}:
        return cast(SubmissionRoute, value)
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
    submission_route: SubmissionRoute
    result_read_source: ResultReadSource
    legacy_worker_endpoints_enabled: bool


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
    submission_route = _parse_submission_route(
        cast(bytes | str | None, redis.get(_SUBMISSION_ROUTE_KEY)),
        cast(SubmissionRoute, settings.cutover_submission_route),
    )
    result_read_source = _parse_result_read_source(
        cast(bytes | str | None, redis.get(_RESULT_READ_SOURCE_KEY)),
        cast(ResultReadSource, settings.cutover_result_read_source),
    )
    legacy_worker_endpoints_enabled = _parse_bool(
        cast(bytes | str | None, redis.get(_LEGACY_WORKER_ENDPOINTS_KEY)),
        settings.legacy_worker_endpoints_enabled,
    )
    return OpsFlags(
        submission_enabled=submission_enabled,
        dequeue_enabled=dequeue_enabled,
        submission_route=submission_route,
        result_read_source=result_read_source,
        legacy_worker_endpoints_enabled=legacy_worker_endpoints_enabled,
    )


def set_ops_flags(
    *,
    submission_enabled: Optional[bool] = None,
    dequeue_enabled: Optional[bool] = None,
    submission_route: Optional[SubmissionRoute] = None,
    result_read_source: Optional[ResultReadSource] = None,
    legacy_worker_endpoints_enabled: Optional[bool] = None,
    redis: Optional[Redis] = None,
) -> OpsFlags:
    redis = redis or get_redis_connection()
    if submission_enabled is not None:
        redis.set(_SUBMISSION_KEY, "1" if submission_enabled else "0")
    if dequeue_enabled is not None:
        redis.set(_DEQUEUE_KEY, "1" if dequeue_enabled else "0")
    if submission_route is not None:
        redis.set(_SUBMISSION_ROUTE_KEY, submission_route)
    if result_read_source is not None:
        redis.set(_RESULT_READ_SOURCE_KEY, result_read_source)
    if legacy_worker_endpoints_enabled is not None:
        redis.set(
            _LEGACY_WORKER_ENDPOINTS_KEY,
            "1" if legacy_worker_endpoints_enabled else "0",
        )
    return get_ops_flags(redis=redis)
