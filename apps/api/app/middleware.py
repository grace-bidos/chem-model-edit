from __future__ import annotations

import logging
import time
from typing import Awaitable, Callable

from fastapi import Request
from starlette.responses import Response

from services.zpe.structured_log import log_event, new_request_id

logger = logging.getLogger(__name__)


def get_request_id(request: Request) -> str:
    request_id = getattr(request.state, "request_id", None)
    return request_id or new_request_id()


def get_tenant_id(request: Request) -> str | None:
    tenant_id = getattr(request.state, "tenant_id", None)
    if isinstance(tenant_id, str):
        value = tenant_id.strip()
        return value or None
    raw = request.headers.get("x-tenant-id")
    if raw is None:
        return None
    value = raw.strip()
    return value or None


async def add_request_context(
    request: Request,
    call_next: Callable[[Request], Awaitable[Response]],
) -> Response:
    """リクエストの相関IDとZPE向けアクセスログを付与する．

    Parameters:
        request: 受信したHTTPリクエスト．
        call_next: 次のミドルウェアまたはエンドポイントを呼ぶコールバック．

    Returns:
        後続処理のレスポンス．
    """
    request_id = request.headers.get("x-request-id") or new_request_id()
    request.state.request_id = request_id
    request.state.tenant_id = get_tenant_id(request)
    start = time.monotonic()
    response = await call_next(request)
    response.headers["x-request-id"] = request_id
    path = request.url.path
    if path.startswith("/api/zpe/"):
        duration_ms = int((time.monotonic() - start) * 1000)
        log_event(
            logger,
            event="zpe_http_request",
            service="control-plane",
            request_id=request_id,
            tenant_id=get_tenant_id(request),
            method=request.method,
            path=path,
            status=response.status_code,
            duration_ms=duration_ms,
        )
    return response
