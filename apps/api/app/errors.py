from __future__ import annotations

import logging
from typing import Any, Optional, cast

from fastapi import HTTPException, Request
from fastapi.exceptions import RequestValidationError
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
from redis.exceptions import RedisError

from app.middleware import get_request_id
from app.schemas.errors import ErrorResponse, ErrorInfo
from services.zpe.structured_log import log_event


logger = logging.getLogger(__name__)

_STATUS_CODE_MAP = {
    400: "bad_request",
    401: "unauthorized",
    403: "forbidden",
    404: "not_found",
    409: "conflict",
    422: "validation_error",
    429: "rate_limited",
    500: "internal_error",
    503: "service_unavailable",
}


def _error_payload(
    *, status_code: int, message: str, details: Optional[Any] = None
) -> dict[str, Any]:
    code = _STATUS_CODE_MAP.get(status_code, "error")
    payload = ErrorResponse(
        error=ErrorInfo(code=code, message=message, details=details)
    )
    return payload.model_dump(by_alias=True, exclude_none=True)


def _sanitize_error_details(details: Any) -> Any:
    if isinstance(details, BaseException):
        return str(details)
    if isinstance(details, dict):
        return {key: _sanitize_error_details(value) for key, value in details.items()}
    if isinstance(details, list):
        return [_sanitize_error_details(value) for value in details]
    return details


def http_exception_handler(_: Request, exc: Exception) -> JSONResponse:
    http_exc = cast(HTTPException, exc)
    message = http_exc.detail
    details: Any = None
    payload = _error_payload(
        status_code=http_exc.status_code,
        message=message,
        details=details,
    )
    return JSONResponse(status_code=http_exc.status_code, content=payload)


def validation_exception_handler(_: Request, exc: Exception) -> JSONResponse:
    validation_exc = cast(RequestValidationError, exc)
    details = jsonable_encoder(_sanitize_error_details(validation_exc.errors()))
    payload = _error_payload(
        status_code=422,
        message="validation error",
        details=details,
    )
    return JSONResponse(status_code=422, content=payload)


def value_error_handler(_: Request, exc: Exception) -> JSONResponse:
    value_exc = cast(ValueError, exc)
    payload = _error_payload(status_code=400, message=str(value_exc))
    return JSONResponse(status_code=400, content=payload)


def overflow_error_handler(_: Request, exc: Exception) -> JSONResponse:
    overflow_exc = cast(OverflowError, exc)
    payload = _error_payload(status_code=400, message=str(overflow_exc))
    return JSONResponse(status_code=400, content=payload)


def redis_error_handler(request: Request, exc: Exception) -> JSONResponse:
    redis_exc = cast(RedisError, exc)
    logger.error("Redis error", exc_info=redis_exc)
    log_event(
        logger,
        event="zpe_redis_error",
        service="control-plane",
        stage="redis",
        status="error",
        request_id=get_request_id(request),
        error_message=str(redis_exc),
    )
    payload = _error_payload(status_code=503, message="redis unavailable")
    return JSONResponse(status_code=503, content=payload)
