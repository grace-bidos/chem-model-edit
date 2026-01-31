from __future__ import annotations

from fastapi import APIRouter, HTTPException, status

from app.schemas.lattice import LatticeConvertRequest, LatticeConvertResponse
from services.lattice import params_to_vectors, vectors_to_params

router = APIRouter(prefix="/api/lattices", tags=["lattices"])


@router.post("/convert", response_model=LatticeConvertResponse)
async def convert_lattice(
    request: LatticeConvertRequest,
) -> LatticeConvertResponse:
    if request.from_ == "vectors":
        if request.lattice is None:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
                detail="lattice is required when from=vectors",
            )
        params = vectors_to_params(request.lattice)
        return LatticeConvertResponse(
            lattice=request.lattice, params=params, unit=request.unit
        )
    if request.params is None:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail="params is required when from=params",
        )
    lattice = params_to_vectors(request.params)
    return LatticeConvertResponse(lattice=lattice, params=request.params, unit=request.unit)
