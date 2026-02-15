from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Response, status

from services.aiida_runtime_checks import probe_managed_aiida_runtime

router = APIRouter(prefix="/api", tags=["health"])


@router.get("/health")
async def health() -> dict[str, Any]:
    aiida_check = probe_managed_aiida_runtime(check_kind="health")
    overall = "ok" if aiida_check.status in {"ok", "skipped"} else "degraded"
    return {
        "status": overall,
        "checks": {
            "managed_aiida_runtime": aiida_check.as_dict(),
        },
    }


@router.get("/ready")
async def ready(response: Response) -> dict[str, Any]:
    aiida_check = probe_managed_aiida_runtime(check_kind="ready")
    is_ready = aiida_check.status in {"ok", "skipped"}
    if not is_ready:
        response.status_code = status.HTTP_503_SERVICE_UNAVAILABLE
    return {
        "status": "ready" if is_ready else "not_ready",
        "checks": {
            "managed_aiida_runtime": aiida_check.as_dict(),
        },
    }
