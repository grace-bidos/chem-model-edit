from __future__ import annotations

from dataclasses import dataclass
import os
import time
from typing import Any, Literal
from urllib import error as urlerror
from urllib import request as urlrequest


ManagedAiidaStatus = Literal["ok", "failed", "skipped"]

_ENABLED_KEY = "AIIA_MANAGED_CHECKS_ENABLED"
_BASE_URL_KEY = "AIIA_MANAGED_BASE_URL"
_HEALTH_PATH_KEY = "AIIA_MANAGED_HEALTH_PATH"
_READY_PATH_KEY = "AIIA_MANAGED_READY_PATH"
_TIMEOUT_SECONDS_KEY = "AIIA_MANAGED_TIMEOUT_SECONDS"
_TOKEN_KEY = "AIIA_MANAGED_BEARER_TOKEN"


def _truthy(value: str | None) -> bool:
    if value is None:
        return False
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _timeout_seconds() -> int:
    raw = os.getenv(_TIMEOUT_SECONDS_KEY, "3").strip()
    try:
        parsed = int(raw)
    except ValueError:
        return 3
    return parsed if parsed > 0 else 3


@dataclass(frozen=True)
class ManagedAiidaCheck:
    status: ManagedAiidaStatus
    detail: str
    url: str | None = None
    error: str | None = None
    http_status: int | None = None
    latency_ms: int | None = None

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "status": self.status,
            "detail": self.detail,
        }
        if self.url is not None:
            payload["url"] = self.url
        if self.error is not None:
            payload["error"] = self.error
        if self.http_status is not None:
            payload["http_status"] = self.http_status
        if self.latency_ms is not None:
            payload["latency_ms"] = self.latency_ms
        return payload


def probe_managed_aiida_runtime(*, check_kind: Literal["health", "ready"]) -> ManagedAiidaCheck:
    """Probe managed AiiDA runtime according to the selected check kind."""

    enabled = _truthy(os.getenv(_ENABLED_KEY))
    if not enabled:
        return ManagedAiidaCheck(
            status="skipped",
            detail="managed AiiDA checks disabled",
        )

    base_url = (os.getenv(_BASE_URL_KEY) or "").strip().rstrip("/")
    if not base_url:
        return ManagedAiidaCheck(
            status="failed",
            detail=f"{_BASE_URL_KEY} is required when {_ENABLED_KEY}=true",
            error="misconfigured",
        )

    if check_kind == "health":
        path = os.getenv(_HEALTH_PATH_KEY, "/health")
    else:
        path = os.getenv(_READY_PATH_KEY, "/ready")
    target_url = f"{base_url}/{path.lstrip('/')}"

    headers = {"Accept": "application/json"}
    token = (os.getenv(_TOKEN_KEY) or "").strip()
    if token:
        headers["Authorization"] = f"Bearer {token}"
    request = urlrequest.Request(target_url, headers=headers, method="GET")

    timeout = _timeout_seconds()
    started = time.perf_counter()
    try:
        with urlrequest.urlopen(request, timeout=timeout) as response:
            http_status = int(response.status)
            latency_ms = int((time.perf_counter() - started) * 1000)
            if 200 <= http_status < 300:
                return ManagedAiidaCheck(
                    status="ok",
                    detail="managed AiiDA runtime reachable",
                    url=target_url,
                    http_status=http_status,
                    latency_ms=latency_ms,
                )
            return ManagedAiidaCheck(
                status="failed",
                detail=f"managed AiiDA returned HTTP {http_status}",
                url=target_url,
                error="upstream_unhealthy",
                http_status=http_status,
                latency_ms=latency_ms,
            )
    except urlerror.HTTPError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return ManagedAiidaCheck(
            status="failed",
            detail=f"managed AiiDA returned HTTP {exc.code}",
            url=target_url,
            error="upstream_unhealthy",
            http_status=int(exc.code),
            latency_ms=latency_ms,
        )
    except urlerror.URLError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return ManagedAiidaCheck(
            status="failed",
            detail=f"managed AiiDA endpoint unreachable: {exc.reason}",
            url=target_url,
            error="unreachable",
            latency_ms=latency_ms,
        )
