from __future__ import annotations

import logging
from typing import Any

from fastapi import APIRouter, HTTPException

from app.schemas.supercell import (
    SupercellBuildMeta,
    SupercellBuildRequest,
    SupercellBuildResponse,
)
from services.parse import structure_from_ase
from services.structures import get_structure_entry, register_structure_atoms
from services.supercell import build_supercell_from_grid

router = APIRouter(prefix="/api/supercells", tags=["supercells"])
logger = logging.getLogger(__name__)


def _cell_has_lattice(cell: Any) -> bool:
    try:
        lengths = cell.lengths()
    except AttributeError:
        return False
    return any(length > 1.0e-8 for length in lengths)


def _cells_match(cell_a: Any, cell_b: Any, *, tol: float = 1.0e-6) -> bool:
    try:
        arr_a = cell_a.array
        arr_b = cell_b.array
    except AttributeError:
        return False
    for row_a, row_b in zip(arr_a, arr_b, strict=True):
        for value_a, value_b in zip(row_a, row_b, strict=True):
            if abs(float(value_a) - float(value_b)) > tol:
                return False
    return True


@router.post("/builds", response_model=SupercellBuildResponse)
async def supercell_build(request: SupercellBuildRequest) -> SupercellBuildResponse:
    try:
        base_entry = get_structure_entry(request.base_structure_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc

    base_atoms = base_entry.atoms
    base_cell = base_atoms.get_cell()
    if not _cell_has_lattice(base_cell):
        raise HTTPException(status_code=400, detail="base structure has no lattice")

    tiles: dict[str, Any] = {}
    structure_ids_used: list[str] = []
    for row in request.grid.tiles:
        for structure_id in row:
            if structure_id in tiles:
                continue
            structure_ids_used.append(structure_id)
            try:
                entry = get_structure_entry(structure_id)
            except KeyError as exc:
                raise HTTPException(
                    status_code=404, detail="Structure not found"
                ) from exc
            tiles[structure_id] = entry.atoms

    options = request.options
    validate_lattice = options.validate_lattice if options else "none"
    if validate_lattice != "none":
        for structure_id, atoms in tiles.items():
            if structure_id == request.base_structure_id:
                continue
            cell = atoms.get_cell()
            lattice_ok = _cell_has_lattice(cell) and _cells_match(base_cell, cell)
            if lattice_ok:
                continue
            if validate_lattice == "error":
                raise HTTPException(
                    status_code=400,
                    detail=f"lattice mismatch for structure {structure_id}",
                )
            logger.warning(
                "supercell.build lattice mismatch: %s vs base %s",
                structure_id,
                request.base_structure_id,
            )

    check_overlap = bool(options.check_overlap) if options else False
    overlap_tolerance = options.overlap_tolerance if options else None
    atoms_out, overlap_count = build_supercell_from_grid(
        base_atoms,
        request.grid,
        tiles,
        check_overlap=check_overlap,
        overlap_tolerance=overlap_tolerance,
    )

    structure_id = register_structure_atoms(atoms_out, source="supercell-build")
    include_structure = bool(request.output.include_structure) if request.output else False
    structure = structure_from_ase(atoms_out) if include_structure else None
    meta = SupercellBuildMeta(
        rows=request.grid.rows,
        cols=request.grid.cols,
        tile_count=request.grid.rows * request.grid.cols,
        overlap_count=overlap_count if check_overlap else None,
        base_structure_id=request.base_structure_id,
        structure_ids_used=structure_ids_used,
    )
    return SupercellBuildResponse(
        structure_id=structure_id,
        structure=structure,
        meta=meta,
    )
