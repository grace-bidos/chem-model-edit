from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Response

from app.schemas.structures import (
    StructureCreateRequest,
    StructureCreateResponse,
    StructureExportRequest,
    StructureExportResponse,
    StructureGetResponse,
    StructureParseRequest,
    StructureParseResponse,
)
from services.export import export_qe_in
from services.parse import parse_qe_in, structure_from_ase
from services.structures import (
    create_structure_from_qe,
    get_structure_cif,
    get_structure_entry,
)

router = APIRouter(prefix="/api/structures", tags=["structures"])


@router.post("/parse", response_model=StructureParseResponse)
async def parse_structure(request: StructureParseRequest) -> StructureParseResponse:
    structure = parse_qe_in(request.content)
    return StructureParseResponse(structure=structure)


@router.post("", response_model=StructureCreateResponse)
async def create_structure(request: StructureCreateRequest) -> StructureCreateResponse:
    structure_id, structure, source, params = create_structure_from_qe(
        request.content
    )
    return StructureCreateResponse(
        id=structure_id,
        structure=structure,
        source=source,
        params=params,
        raw_input=request.content,
    )


@router.get("/{structure_id}", response_model=StructureGetResponse)
async def get_structure(structure_id: str) -> StructureGetResponse:
    try:
        entry = get_structure_entry(structure_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc
    structure = structure_from_ase(entry.atoms)
    return StructureGetResponse(
        structure=structure,
        params=entry.params,
        raw_input=entry.raw_input,
        source=entry.source,
    )


@router.get("/{structure_id}/view")
async def view_structure(
    structure_id: str,
    format: str = Query("cif"),
) -> Response:
    if format != "cif":
        raise HTTPException(status_code=400, detail="Unsupported format")
    try:
        cif = get_structure_cif(structure_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc
    return Response(content=cif, media_type="chemical/x-cif")


@router.post("/export", response_model=StructureExportResponse)
async def export_structure(
    request: StructureExportRequest,
) -> StructureExportResponse:
    content = export_qe_in(request.structure)
    return StructureExportResponse(content=content)
