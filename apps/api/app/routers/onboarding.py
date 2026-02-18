from __future__ import annotations

from fastapi import APIRouter, Request

from app.deps import require_admin
from app.schemas.errors import ErrorResponse
from app.schemas.onboarding import OnboardingDryRunRequest, OnboardingDryRunResponse
from services.onboarding_manifest import run_onboarding_dry_run

router = APIRouter(prefix="/api/onboarding", tags=["onboarding"])


@router.post(
    "/dry-run",
    response_model=OnboardingDryRunResponse,
    responses={
        401: {
            "model": ErrorResponse,
            "description": "Unauthorized (admin required)",
        }
    },
)
async def onboarding_dry_run(
    request: OnboardingDryRunRequest,
    raw: Request,
) -> OnboardingDryRunResponse:
    require_admin(raw)
    return run_onboarding_dry_run(request)
