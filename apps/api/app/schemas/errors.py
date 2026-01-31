from __future__ import annotations

from typing import Any, Optional

from .base import ApiModel


class ErrorInfo(ApiModel):
    code: str
    message: str
    details: Optional[Any] = None


class ErrorResponse(ApiModel):
    error: ErrorInfo
