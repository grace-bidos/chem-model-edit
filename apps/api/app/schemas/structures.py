from __future__ import annotations

from typing import Optional

from .base import ApiModel
from .common import QeParameters, Structure


class StructureParseRequest(ApiModel):
    content: str
    format: Optional[str] = None


class StructureParseResponse(ApiModel):
    structure: Structure


class StructureCreateRequest(ApiModel):
    content: str
    format: Optional[str] = None


class StructureCreateResponse(ApiModel):
    id: str
    structure: Structure
    source: str
    params: Optional[QeParameters] = None
    raw_input: Optional[str] = None


class StructureGetResponse(ApiModel):
    structure: Structure
    params: Optional[QeParameters] = None
    raw_input: Optional[str] = None
    source: Optional[str] = None


class StructureExportRequest(ApiModel):
    structure: Structure
    format: Optional[str] = None


class StructureExportResponse(ApiModel):
    content: str
