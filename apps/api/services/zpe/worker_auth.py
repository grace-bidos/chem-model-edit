from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import hashlib
import secrets
from typing import Any, Optional, cast

from redis import Redis

from .queue import get_redis_connection
from .settings import get_zpe_settings


_TOKEN_PREFIX = "zpe:token:"
_REVOKED_SET = "zpe:revoked_tokens"
_WORKER_INDEX_PREFIX = "zpe:worker_tokens:"


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


def _now_ts() -> int:
    return int(_now().timestamp())


def _hash_token(token: str) -> str:
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


@dataclass
class WorkerToken:
    token: str
    token_hash: str
    worker_id: str
    tenant_id: Optional[str]
    workspace_id: Optional[str]
    expires_at: str
    ttl_seconds: int


@dataclass
class WorkerPrincipal:
    worker_id: str
    tenant_id: Optional[str] = None
    workspace_id: Optional[str] = None


def _normalize_scope(value: bytes | str | None) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, bytes):
        raw = value.decode("utf-8")
    else:
        raw = value
    normalized = raw.strip()
    return normalized or None


class WorkerTokenStore:
    def __init__(self, redis: Optional[Redis] = None) -> None:
        self.redis = redis or get_redis_connection()

    def create_token(
        self,
        worker_id: str,
        *,
        label: Optional[str] = None,
        tenant_id: Optional[str] = None,
        workspace_id: Optional[str] = None,
    ) -> WorkerToken:
        settings = get_zpe_settings()
        ttl = int(settings.worker_token_ttl_seconds)
        if ttl <= 0:
            raise ValueError("worker_token_ttl_seconds must be >= 1")
        token = secrets.token_urlsafe(32)
        token_hash = _hash_token(token)
        expires_at = (_now() + timedelta(seconds=ttl)).isoformat()
        key = f"{_TOKEN_PREFIX}{token_hash}"
        payload = {
            "worker_id": worker_id,
            "tenant_id": tenant_id or "",
            "workspace_id": workspace_id or "",
            "created_at": _now_iso(),
            "expires_at": expires_at,
            "label": label or "",
            "revoked_at": "",
        }
        redis_any = cast(Any, self.redis)
        pipe = redis_any.pipeline(transaction=True)
        pipe.hset(key, mapping=payload)
        pipe.expire(key, ttl)
        pipe.sadd(f"{_WORKER_INDEX_PREFIX}{worker_id}", token_hash)
        pipe.expire(f"{_WORKER_INDEX_PREFIX}{worker_id}", ttl)
        pipe.execute()
        return WorkerToken(
            token=token,
            token_hash=token_hash,
            worker_id=worker_id,
            tenant_id=tenant_id,
            workspace_id=workspace_id,
            expires_at=expires_at,
            ttl_seconds=ttl,
        )

    def validate(self, token: str) -> WorkerPrincipal:
        token_hash = _hash_token(token)
        self._cleanup_revoked()
        if self.redis.zscore(_REVOKED_SET, token_hash) is not None:
            raise PermissionError("token revoked")
        key = f"{_TOKEN_PREFIX}{token_hash}"
        redis_any = cast(Any, self.redis)
        data = cast(dict[bytes, bytes], redis_any.hgetall(key))
        if not data:
            raise PermissionError("token invalid or expired")
        worker_id = data.get(b"worker_id", b"").decode("utf-8")
        if not worker_id:
            raise PermissionError("token invalid")
        return WorkerPrincipal(
            worker_id=worker_id,
            tenant_id=_normalize_scope(data.get(b"tenant_id")),
            workspace_id=_normalize_scope(data.get(b"workspace_id")),
        )

    def revoke_tokens_for_worker(self, worker_id: str) -> int:
        settings = get_zpe_settings()
        index_key = f"{_WORKER_INDEX_PREFIX}{worker_id}"
        redis_any = cast(Any, self.redis)
        token_hashes = cast(set[bytes], redis_any.smembers(index_key))
        if not token_hashes:
            return 0
        pipe = redis_any.pipeline(transaction=True)
        revoked_at = _now_ts()
        for token_hash_raw in token_hashes:
            token_hash = token_hash_raw.decode("utf-8")
            pipe.zadd(_REVOKED_SET, {token_hash: revoked_at})
            pipe.hset(
                f"{_TOKEN_PREFIX}{token_hash}", mapping={"revoked_at": _now_iso()}
            )
        pipe.expire(_REVOKED_SET, int(settings.worker_token_ttl_seconds * 2))
        pipe.execute()
        return len(token_hashes)

    def _cleanup_revoked(self) -> None:
        ttl = int(get_zpe_settings().worker_token_ttl_seconds)
        if ttl <= 0:
            return
        cutoff = _now_ts() - ttl
        self.redis.zremrangebyscore(_REVOKED_SET, 0, cutoff)


def get_worker_token_store() -> WorkerTokenStore:
    return WorkerTokenStore()
