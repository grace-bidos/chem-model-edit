from __future__ import annotations

from fastapi import APIRouter

from app.schemas.transforms import DeltaTransplantRequest, DeltaTransplantResponse
from services.transplant import transplant_delta

router = APIRouter(prefix="/api/transforms", tags=["transforms"])


@router.post("/delta-transplant", response_model=DeltaTransplantResponse)
async def delta_transplant(
    request: DeltaTransplantRequest,
) -> DeltaTransplantResponse:
    content = transplant_delta(
        request.small_in,
        request.small_out,
        request.large_in,
    )
    return DeltaTransplantResponse(content=content)
