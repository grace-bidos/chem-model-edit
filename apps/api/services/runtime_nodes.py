from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import json
import os
import secrets
from typing import Any, Optional, cast
from urllib import error as urlerror
from urllib import request as urlrequest
from uuid import uuid4


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _now_iso() -> str:
    return _now().isoformat()


def _expires_iso(ttl_seconds: int) -> str:
    return (_now() + timedelta(seconds=ttl_seconds)).isoformat()


def _normalize_queue_name(queue_name: str) -> str:
    normalized = queue_name.strip()
    if not normalized:
        raise ValueError("queue_name must be a non-empty string")
    return normalized


@dataclass(frozen=True)
class RuntimeTarget:
    target_id: str
    tenant_id: str
    user_id: str
    queue_name: str
    server_id: str
    registered_at: str
    name: Optional[str] = None


@dataclass(frozen=True)
class RuntimeJoinToken:
    token: str
    tenant_id: str
    owner_user_id: str
    queue_name: str
    created_at: str
    expires_at: str
    ttl_seconds: int
    node_name_hint: Optional[str] = None
    label: Optional[str] = None


@dataclass(frozen=True)
class RuntimeNodeRegistration:
    server_id: str
    target_id: str
    queue_name: str
    registered_at: str
    name: Optional[str]


class RuntimeNodeStore:
    def __init__(self) -> None:
        convex_url = os.getenv("CONVEX_URL", "").strip()
        deploy_key = os.getenv("CONVEX_DEPLOY_KEY", "").strip()
        if not convex_url or not deploy_key:
            raise RuntimeError("CONVEX_URL and CONVEX_DEPLOY_KEY are required")
        self._convex_url = convex_url.rstrip("/")
        self._deploy_key = deploy_key
        timeout_raw = os.getenv("RUNTIME_NODE_STORE_TIMEOUT_SECONDS", "10")
        try:
            timeout = int(timeout_raw)
        except ValueError:
            timeout = 10
        self._timeout_seconds = max(1, timeout)

    def _request(self, method: str, path: str, args: dict[str, Any]) -> Any:
        body = json.dumps(
            {"path": path, "args": args},
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("utf-8")
        req = urlrequest.Request(
            f"{self._convex_url}/api/{method}",
            method="POST",
            data=body,
            headers={
                "Authorization": f"Convex {self._deploy_key}",
                "Content-Type": "application/json",
            },
        )
        try:
            with urlrequest.urlopen(req, timeout=self._timeout_seconds) as resp:
                raw = resp.read().decode("utf-8", errors="replace")
        except urlerror.HTTPError as exc:
            detail = exc.read().decode("utf-8", errors="replace")
            raise RuntimeError(f"convex {method} failed: HTTP {exc.code} {detail}") from exc
        except urlerror.URLError as exc:
            raise RuntimeError(f"convex {method} failed") from exc

        payload_obj: object = json.loads(raw) if raw.strip() else {}
        if not isinstance(payload_obj, dict):
            return {}
        payload = cast(dict[str, Any], payload_obj)
        return payload.get("value", payload)

    def create_join_token(
        self,
        *,
        tenant_id: str,
        owner_user_id: str,
        queue_name: str,
        ttl_seconds: int,
        node_name_hint: str | None = None,
    ) -> RuntimeJoinToken:
        if ttl_seconds <= 0:
            raise ValueError("ttl_seconds must be >= 1")
        token = secrets.token_urlsafe(32)
        now = _now_iso()
        expires_at = _expires_iso(ttl_seconds)
        normalized_queue = _normalize_queue_name(queue_name)
        self._request(
            "mutation",
            "runtimeTargets:createJoinToken",
            {
                "token": token,
                "tenant_id": tenant_id,
                "owner_user_id": owner_user_id,
                "queue_name": normalized_queue,
                "created_at": now,
                "expires_at": expires_at,
                "node_name_hint": node_name_hint,
                "label": normalized_queue,
            },
        )
        return RuntimeJoinToken(
            token=token,
            tenant_id=tenant_id,
            owner_user_id=owner_user_id,
            queue_name=normalized_queue,
            created_at=now,
            expires_at=expires_at,
            ttl_seconds=ttl_seconds,
            node_name_hint=node_name_hint,
            label=normalized_queue,
        )

    def consume_join_token(self, token: str) -> RuntimeJoinToken:
        value = self._request(
            "mutation",
            "runtimeTargets:consumeJoinToken",
            {"token": token, "consumed_at": _now_iso()},
        )
        if not isinstance(value, dict):
            raise KeyError("join token not found")
        value_map = cast(dict[str, Any], value)
        if not value_map.get("ok"):
            raise KeyError(str(value_map.get("reason", "join token not found")))
        join_payload = value_map.get("join_token")
        if not isinstance(join_payload, dict):
            raise RuntimeError("invalid join token payload")
        payload = cast(dict[str, Any], join_payload)
        return RuntimeJoinToken(
            token=cast(str, payload.get("token", token)),
            tenant_id=cast(str, payload["tenant_id"]),
            owner_user_id=cast(str, payload["owner_user_id"]),
            queue_name=cast(str, payload["queue_name"]),
            created_at=cast(str, payload["created_at"]),
            expires_at=cast(str, payload["expires_at"]),
            ttl_seconds=0,
            node_name_hint=cast(Optional[str], payload.get("node_name_hint")),
            label=cast(Optional[str], payload.get("label")),
        )

    def add_target(
        self,
        *,
        tenant_id: str,
        user_id: str,
        queue_name: str,
        server_id: str,
        name: str | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> RuntimeTarget:
        target_id = f"qt-{uuid4().hex}"
        registered_at = _now_iso()
        normalized_queue = _normalize_queue_name(queue_name)
        self._request(
            "mutation",
            "runtimeTargets:addRuntimeTarget",
            {
                "tenant_id": tenant_id,
                "user_id": user_id,
                "target_id": target_id,
                "queue_name": normalized_queue,
                "server_id": server_id,
                "registered_at": registered_at,
                "name": name,
                "metadata": metadata,
            },
        )
        return RuntimeTarget(
            target_id=target_id,
            tenant_id=tenant_id,
            user_id=user_id,
            queue_name=normalized_queue,
            server_id=server_id,
            registered_at=registered_at,
            name=name,
        )

    def list_targets(self, tenant_id: str, user_id: str) -> list[RuntimeTarget]:
        value = self._request(
            "query",
            "runtimeTargets:listTargets",
            {"tenant_id": tenant_id, "user_id": user_id},
        )
        if not isinstance(value, list):
            return []
        targets: list[RuntimeTarget] = []
        for item in value:
            if not isinstance(item, dict):
                continue
            payload = cast(dict[str, Any], item)
            target_id = payload.get("target_id")
            queue_name = payload.get("queue_name")
            server_id = payload.get("server_id")
            registered_at = payload.get("registered_at")
            row_user_id = payload.get("user_id")
            row_tenant_id = payload.get("tenant_id")
            if (
                isinstance(target_id, str)
                and isinstance(queue_name, str)
                and isinstance(server_id, str)
                and isinstance(registered_at, str)
                and isinstance(row_user_id, str)
                and isinstance(row_tenant_id, str)
            ):
                targets.append(
                    RuntimeTarget(
                        target_id=target_id,
                        tenant_id=row_tenant_id,
                        user_id=row_user_id,
                        queue_name=queue_name,
                        server_id=server_id,
                        registered_at=registered_at,
                        name=cast(Optional[str], payload.get("name")),
                    )
                )
        return targets

    def get_target_by_id(self, target_id: str) -> Optional[RuntimeTarget]:
        value = self._request(
            "query",
            "runtimeTargets:getTargetById",
            {"target_id": target_id},
        )
        if not isinstance(value, dict):
            return None
        payload = cast(dict[str, Any], value)
        queue_name = payload.get("queue_name")
        server_id = payload.get("server_id")
        registered_at = payload.get("registered_at")
        user_id = payload.get("user_id")
        tenant_id = payload.get("tenant_id")
        if (
            not isinstance(queue_name, str)
            or not isinstance(server_id, str)
            or not isinstance(registered_at, str)
            or not isinstance(user_id, str)
            or not isinstance(tenant_id, str)
        ):
            return None
        return RuntimeTarget(
            target_id=target_id,
            tenant_id=tenant_id,
            user_id=user_id,
            queue_name=queue_name,
            server_id=server_id,
            registered_at=registered_at,
            name=cast(Optional[str], payload.get("name")),
        )

    def get_active_target(self, tenant_id: str, user_id: str) -> Optional[RuntimeTarget]:
        value = self._request(
            "query",
            "runtimeTargets:getActiveTarget",
            {"tenant_id": tenant_id, "user_id": user_id},
        )
        if not isinstance(value, dict):
            return None
        payload = cast(dict[str, Any], value)
        target_id = payload.get("target_id")
        queue_name = payload.get("queue_name")
        server_id = payload.get("server_id")
        registered_at = payload.get("registered_at")
        row_user_id = payload.get("user_id")
        row_tenant_id = payload.get("tenant_id")
        if (
            not isinstance(target_id, str)
            or not isinstance(queue_name, str)
            or not isinstance(server_id, str)
            or not isinstance(registered_at, str)
            or not isinstance(row_user_id, str)
            or not isinstance(row_tenant_id, str)
        ):
            return None
        return RuntimeTarget(
            target_id=target_id,
            tenant_id=row_tenant_id,
            user_id=row_user_id,
            queue_name=queue_name,
            server_id=server_id,
            registered_at=registered_at,
            name=cast(Optional[str], payload.get("name")),
        )

    def set_active_target(self, tenant_id: str, user_id: str, target_id: str) -> None:
        self._request(
            "mutation",
            "runtimeTargets:setActiveTarget",
            {
                "tenant_id": tenant_id,
                "user_id": user_id,
                "target_id": target_id,
                "updated_at": _now_iso(),
            },
        )

    def ensure_target_owner(self, tenant_id: str, user_id: str, target_id: str) -> RuntimeTarget:
        target = self.get_target_by_id(target_id)
        if not target or target.user_id != user_id or target.tenant_id != tenant_id:
            raise KeyError("target not found")
        return target


_runtime_node_store: RuntimeNodeStore | None = None


def get_runtime_node_store() -> RuntimeNodeStore:
    global _runtime_node_store
    if _runtime_node_store is None:
        _runtime_node_store = RuntimeNodeStore()
    return _runtime_node_store
