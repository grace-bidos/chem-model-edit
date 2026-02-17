from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Optional, cast

from redis import Redis

from .job_meta import get_job_meta_store
from .queue import get_redis_connection
from .settings import get_zpe_settings


_OWNER_PREFIX = "zpe:job:owner:"


@dataclass(frozen=True)
class JobOwnerRecord:
    user_id: str
    tenant_id: str | None


class JobOwnerStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def set_owner(self, job_id: str, user_id: str, tenant_id: str) -> None:
        ttl = get_zpe_settings().result_ttl_seconds
        payload = json.dumps(
            {"user_id": user_id, "tenant_id": tenant_id},
            ensure_ascii=True,
            separators=(",", ":"),
        )
        ok = self.redis.setex(f"{_OWNER_PREFIX}{job_id}", ttl, payload)
        if not ok:
            raise RuntimeError("failed to set job owner")

    def get_owner_record(self, job_id: str) -> Optional[JobOwnerRecord]:
        # Runtime ownership should come from execution metadata first.
        meta = get_job_meta_store(redis=self.redis).get_meta(job_id)
        meta_user_id = meta.get("user_id")
        if isinstance(meta_user_id, str) and meta_user_id:
            meta_tenant_id = meta.get("tenant_id")
            if meta_tenant_id is not None and not isinstance(meta_tenant_id, str):
                return None
            return JobOwnerRecord(
                user_id=meta_user_id,
                tenant_id=meta_tenant_id,
            )

        # Backward-compatible fallback for pre-cutover legacy keys.
        raw = cast(Optional[bytes], self.redis.get(f"{_OWNER_PREFIX}{job_id}"))
        if not raw:
            return None
        decoded = raw.decode("utf-8")
        if decoded.startswith("{"):
            try:
                payload = cast(dict[str, object], json.loads(decoded))
            except json.JSONDecodeError:
                return None
            user_id = payload.get("user_id")
            tenant_id = payload.get("tenant_id")
            if not isinstance(user_id, str) or not user_id:
                return None
            if tenant_id is not None and not isinstance(tenant_id, str):
                return None
            return JobOwnerRecord(user_id=user_id, tenant_id=tenant_id)
        return JobOwnerRecord(user_id=decoded, tenant_id=None)

    def get_owner(self, job_id: str) -> Optional[str]:
        owner = self.get_owner_record(job_id)
        return owner.user_id if owner else None


def get_job_owner_store() -> JobOwnerStore:
    return JobOwnerStore()
