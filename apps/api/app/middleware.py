from __future__ import annotations

import logging
import time

from fastapi import Request

from services.zpe.structured_log import log_event, new_request_id

logger = logging.getLogger(__name__)


def get_request_id(request: Request) -> str:
    request_id = getattr(request.state, "request_id", None)
    return request_id or new_request_id()


async def add_request_context(request: Request, call_next):
    request_id = request.headers.get("x-request-id") or new_request_id()
    request.state.request_id = request_id
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
            method=request.method,
            path=path,
            status=response.status_code,
            duration_ms=duration_ms,
        )
    return response
