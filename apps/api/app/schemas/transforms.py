from __future__ import annotations

from .base import ApiModel


class DeltaTransplantRequest(ApiModel):
    small_in: str
    small_out: str
    large_in: str


class DeltaTransplantResponse(ApiModel):
    content: str
