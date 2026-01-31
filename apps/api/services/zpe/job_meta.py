from __future__ import annotations

import json
from typing import Any, Dict, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_META_PREFIX = "zpe:job:meta:"


class JobMetaStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def set_meta(self, job_id: str, meta: Dict[str, Any]) -> None:
        settings = get_zpe_settings()
        key = f"{_META_PREFIX}{job_id}"
        payload = json.dumps(meta, ensure_ascii=True)
        ok = self.redis.setex(key, settings.result_ttl_seconds, payload)
        if not ok:
            raise RuntimeError("failed to set job meta")

    def get_meta(self, job_id: str) -> Dict[str, Any]:
        raw = cast(Optional[bytes], self.redis.get(f"{_META_PREFIX}{job_id}"))
        if not raw:
            return {}
        try:
            return cast(Dict[str, Any], json.loads(raw))
        except json.JSONDecodeError:
            return {}


def get_job_meta_store() -> JobMetaStore:
    return JobMetaStore()
