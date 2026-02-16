from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import time
from typing import Any, Literal, Optional, cast
from uuid import uuid4

from redis import Redis
from redis.exceptions import WatchError

from .queue import get_redis_connection
from .settings import get_zpe_settings

_SUBMIT_PREFIX = "zpe:submit:idempotency:"


@dataclass(frozen=True)
class SubmitIdempotencyRecord:
    job_id: str
    request_fingerprint: str
    requested_queue_name: str
    resolved_queue_name: str


@dataclass(frozen=True)
class SubmitIdempotencyClaim:
    state: Literal["claimed", "ready", "pending"]
    claim_token: str | None = None
    record: SubmitIdempotencyRecord | None = None


class SubmitIdempotencyStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def _key(self, *, user_id: str, tenant_id: str, request_id: str) -> str:
        _ = user_id
        return f"{_SUBMIT_PREFIX}{tenant_id}:{request_id}"

    def _legacy_key(self, *, user_id: str, tenant_id: str, request_id: str) -> str:
        return f"{_SUBMIT_PREFIX}{tenant_id}:{user_id}:{request_id}"

    def _resolve_record_key(self, *, user_id: str, tenant_id: str, request_id: str) -> str | None:
        primary = self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        if self.redis.exists(primary):
            return primary
        legacy = self._legacy_key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        if self.redis.exists(legacy):
            return legacy
        return None

    def get_record(
        self, *, user_id: str, tenant_id: str, request_id: str
    ) -> SubmitIdempotencyRecord | None:
        key = self._resolve_record_key(
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
        )
        if key is None:
            return None
        raw = cast(
            Optional[bytes],
            self.redis.get(key),
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
                "state": "ready",
            },
            ensure_ascii=True,
            separators=(",", ":"),
        )
        ok = self.redis.setex(key, settings.result_ttl_seconds, payload)
        if not ok:
            raise RuntimeError("failed to remember submit idempotency record")

    def claim_or_get(
        self, *, user_id: str, tenant_id: str, request_id: str, request_fingerprint: str
    ) -> SubmitIdempotencyClaim:
        settings = get_zpe_settings()
        claim_token = uuid4().hex
        legacy_key = self._legacy_key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        legacy_current = self._load_raw(legacy_key)
        if legacy_current:
            legacy_state = legacy_current.get("state")
            legacy_fingerprint = legacy_current.get("request_fingerprint")
            if legacy_fingerprint != request_fingerprint:
                raise ValueError("request_id already used for a different submission payload")
            if legacy_state == "ready":
                legacy_record = self.get_record(
                    user_id=user_id,
                    tenant_id=tenant_id,
                    request_id=request_id,
                )
                if legacy_record:
                    return SubmitIdempotencyClaim(state="ready", record=legacy_record)
                raise RuntimeError("submit idempotency record is malformed")
            if legacy_state == "pending":
                return SubmitIdempotencyClaim(state="pending")
        key = self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        pending_payload = json.dumps(
            {
                "state": "pending",
                "request_fingerprint": request_fingerprint,
                "claim_token": claim_token,
            },
            ensure_ascii=True,
            separators=(",", ":"),
        )
        claimed = self.redis.set(
            key,
            pending_payload,
            ex=settings.result_ttl_seconds,
            nx=True,
        )
        if claimed:
            return SubmitIdempotencyClaim(state="claimed", claim_token=claim_token)
        current = self._load_raw(key)
        state = current.get("state")
        existing_fingerprint = current.get("request_fingerprint")
        if existing_fingerprint != request_fingerprint:
            raise ValueError("request_id already used for a different submission payload")
        if state == "ready":
            record = self.get_record(
                user_id=user_id,
                tenant_id=tenant_id,
                request_id=request_id,
            )
            if record:
                return SubmitIdempotencyClaim(state="ready", record=record)
            raise RuntimeError("submit idempotency record is malformed")
        if state == "pending":
            return SubmitIdempotencyClaim(state="pending")
        raise RuntimeError("submit idempotency state is invalid")

    def wait_for_record(
        self,
        *,
        user_id: str,
        tenant_id: str,
        request_id: str,
        request_fingerprint: str,
        timeout_seconds: float = 1.0,
        poll_interval_seconds: float = 0.01,
    ) -> SubmitIdempotencyRecord | None:
        deadline = time.monotonic() + timeout_seconds
        while time.monotonic() < deadline:
            raw = self._load_raw(
                self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
            )
            state = raw.get("state")
            existing_fingerprint = raw.get("request_fingerprint")
            if not state:
                return None
            if existing_fingerprint != request_fingerprint:
                raise ValueError("request_id already used for a different submission payload")
            if state == "ready":
                return self.get_record(
                    user_id=user_id,
                    tenant_id=tenant_id,
                    request_id=request_id,
                )
            time.sleep(poll_interval_seconds)
        return None

    def finalize_claim(
        self,
        *,
        user_id: str,
        tenant_id: str,
        request_id: str,
        claim_token: str,
        record: SubmitIdempotencyRecord,
    ) -> None:
        settings = get_zpe_settings()
        key = self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        ready_payload = json.dumps(
            {
                "state": "ready",
                "job_id": record.job_id,
                "request_fingerprint": record.request_fingerprint,
                "requested_queue_name": record.requested_queue_name,
                "resolved_queue_name": record.resolved_queue_name,
            },
            ensure_ascii=True,
            separators=(",", ":"),
        )
        for _ in range(8):
            with self.redis.pipeline() as pipe:
                try:
                    pipe.watch(key)
                    current = self._load_raw(key, pipe=pipe)
                    if current.get("state") != "pending":
                        pipe.reset()
                        return
                    if current.get("claim_token") != claim_token:
                        pipe.reset()
                        return
                    pipe.multi()
                    pipe.set(key, ready_payload, ex=settings.result_ttl_seconds)
                    pipe.execute()
                    return
                except WatchError:
                    continue
        raise RuntimeError("failed to finalize submit idempotency claim")

    def release_claim(
        self, *, user_id: str, tenant_id: str, request_id: str, claim_token: str
    ) -> None:
        key = self._key(user_id=user_id, tenant_id=tenant_id, request_id=request_id)
        for _ in range(8):
            with self.redis.pipeline() as pipe:
                try:
                    pipe.watch(key)
                    current = self._load_raw(key, pipe=pipe)
                    if current.get("state") != "pending":
                        pipe.reset()
                        return
                    if current.get("claim_token") != claim_token:
                        pipe.reset()
                        return
                    pipe.multi()
                    pipe.delete(key)
                    pipe.execute()
                    return
                except WatchError:
                    continue
        raise RuntimeError("failed to release submit idempotency claim")

    def _load_raw(self, key: str, *, pipe: Any | None = None) -> dict[str, Any]:
        client = pipe if pipe is not None else self.redis
        raw = cast(Optional[bytes], client.get(key))
        if not raw:
            return {}
        try:
            return cast(dict[str, Any], json.loads(raw))
        except json.JSONDecodeError:
            return {}


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
