from __future__ import annotations

from fastapi import APIRouter, HTTPException, Query, Response
from ase import Atoms as ASEAtoms

from app.schemas.structures import (
    StructureCreateRequest,
    StructureCreateResponse,
    StructureExportRequest,
    StructureExportResponse,
    StructureGetResponse,
    StructureParseRequest,
    StructureParseResponse,
)
from app.schemas.common import Structure
from services.cif import atoms_to_cif
from services.export import export_qe_in
from services.parse import parse_qe_in, structure_from_ase
from services.structures import (
    create_structure_from_qe,
    get_structure_cif,
    get_structure_entry,
)

router = APIRouter(prefix="/api/structures", tags=["structures"])


def _ase_from_structure(structure: Structure) -> ASEAtoms:
    if not structure.atoms:
        raise ValueError("原子が空です。")
    symbols = [atom.symbol for atom in structure.atoms]
    positions = [(atom.x, atom.y, atom.z) for atom in structure.atoms]
    if structure.lattice is None:
        return ASEAtoms(symbols=symbols, positions=positions)
    lattice = structure.lattice
    cell = [
        (lattice.a.x, lattice.a.y, lattice.a.z),
        (lattice.b.x, lattice.b.y, lattice.b.z),
        (lattice.c.x, lattice.c.y, lattice.c.z),
    ]
    return ASEAtoms(symbols=symbols, positions=positions, cell=cell, pbc=True)


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
        structure_id=structure_id,
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


@router.post("/cif")
async def export_structure_cif(
    request: StructureExportRequest,
) -> Response:
    try:
        atoms = _ase_from_structure(request.structure)
        cif_text = atoms_to_cif(atoms)
    except Exception as exc:
        raise ValueError("CIFへの変換に失敗しました。") from exc
    return Response(content=cif_text, media_type="chemical/x-cif")
