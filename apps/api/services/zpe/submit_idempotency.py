from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Any, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings

_SUBMIT_PREFIX = "zpe:submit:idempotency:"


@dataclass(frozen=True)
class SubmitIdempotencyRecord:
    job_id: str
    request_fingerprint: str
    requested_queue_name: str
    resolved_queue_name: str


class SubmitIdempotencyStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def _key(self, *, user_id: str, tenant_id: str, request_id: str) -> str:
        return f"{_SUBMIT_PREFIX}{tenant_id}:{user_id}:{request_id}"

    def get_record(
        self, *, user_id: str, tenant_id: str, request_id: str
    ) -> SubmitIdempotencyRecord | None:
        raw = cast(
            Optional[bytes],
            self.redis.get(
                self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
            ),
        )
        if not raw:
            return None
        try:
            payload = cast(dict[str, Any], json.loads(raw))
        except json.JSONDecodeError:
            return None
        job_id = payload.get("job_id")
        request_fingerprint = payload.get("request_fingerprint")
        requested_queue_name = payload.get("requested_queue_name")
        resolved_queue_name = payload.get("resolved_queue_name")
        if (
            not isinstance(job_id, str)
            or not job_id
            or not isinstance(request_fingerprint, str)
            or not request_fingerprint
            or not isinstance(requested_queue_name, str)
            or not requested_queue_name
            or not isinstance(resolved_queue_name, str)
            or not resolved_queue_name
        ):
            return None
        return SubmitIdempotencyRecord(
            job_id=job_id,
            request_fingerprint=request_fingerprint,
            requested_queue_name=requested_queue_name,
            resolved_queue_name=resolved_queue_name,
        )

    def remember(
        self,
        *,
        user_id: str,
        tenant_id: str,
        request_id: str,
        record: SubmitIdempotencyRecord,
    ) -> None:
        settings = get_zpe_settings()
        key = self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        payload = json.dumps(
            {
                "job_id": record.job_id,
                "request_fingerprint": record.request_fingerprint,
                "requested_queue_name": record.requested_queue_name,
                "resolved_queue_name": record.resolved_queue_name,
            },
            ensure_ascii=True,
            separators=(",", ":"),
        )
        ok = self.redis.setex(key, settings.result_ttl_seconds, payload)
        if not ok:
            raise RuntimeError("failed to remember submit idempotency record")


def compute_submit_request_fingerprint(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def get_submit_idempotency_store() -> SubmitIdempotencyStore:
    return SubmitIdempotencyStore()
