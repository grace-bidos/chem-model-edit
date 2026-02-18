from __future__ import annotations

from typing import Any, cast

from redis import Redis

from .settings import get_zpe_settings


def get_redis_connection() -> Redis:
    settings = get_zpe_settings()
    redis_cls = cast(Any, Redis)
    return cast(Redis, redis_cls.from_url(settings.redis_url))


def get_queue(name: str | None = None) -> Any:
    _ = name
    raise RuntimeError("RQ compute queue is retired; use /api/runtime/* integration")


def fetch_job(job_id: str) -> Any:
    _ = job_id
    raise RuntimeError("RQ compute queue is retired; use /api/runtime/* integration")
