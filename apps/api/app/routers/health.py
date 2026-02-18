from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Response, status

from services.aiida_runtime_checks import (
    probe_managed_aiida_runtime,
    probe_slurm_real_adapter_preconditions,
    probe_user_managed_deep_readiness,
)
from services.runtime_settings import get_runtime_settings

router = APIRouter(prefix="/api", tags=["health"])


def _runtime_gateway_contract_check() -> dict[str, Any]:
    settings = get_runtime_settings()
    required = {
        "RUNTIME_COMMAND_SUBMIT_URL": settings.command_submit_url,
        "RUNTIME_COMMAND_EVENT_URL_TEMPLATE": settings.command_event_url_template,
        "RUNTIME_READ_STATUS_URL_TEMPLATE": settings.read_status_url_template,
        "RUNTIME_READ_DETAIL_URL_TEMPLATE": settings.read_detail_url_template,
        "RUNTIME_READ_PROJECTION_URL_TEMPLATE": settings.read_projection_url_template,
    }
    missing = [key for key, value in required.items() if not (value and value.strip())]
    if missing:
        return {
            "status": "failed",
            "detail": "runtime gateway settings are incomplete",
            "missing": missing,
        }
    return {
        "status": "ok",
        "detail": "runtime gateway settings are configured",
    }


@router.get("/health")
def health() -> dict[str, Any]:
    aiida_check = probe_managed_aiida_runtime(check_kind="health")
    deep_readiness_check = probe_user_managed_deep_readiness(check_kind="health")
    real_adapter_check = probe_slurm_real_adapter_preconditions(check_kind="health")
    runtime_gateway_check = _runtime_gateway_contract_check()
    overall = (
        "ok"
        if aiida_check.status in {"ok", "skipped"}
        and deep_readiness_check.status in {"ok", "skipped"}
        and real_adapter_check.status in {"ok", "skipped"}
        and runtime_gateway_check["status"] == "ok"
        else "degraded"
    )
    return {
        "status": overall,
        "checks": {
            "managed_aiida_runtime": aiida_check.as_dict(),
            "user_managed_deep_readiness": deep_readiness_check.as_dict(),
            "slurm_real_adapter_preconditions": real_adapter_check.as_dict(),
            "runtime_gateway_contract": runtime_gateway_check,
        },
    }


@router.get("/ready")
def ready(response: Response) -> dict[str, Any]:
    aiida_check = probe_managed_aiida_runtime(check_kind="ready")
    deep_readiness_check = probe_user_managed_deep_readiness(check_kind="ready")
    real_adapter_check = probe_slurm_real_adapter_preconditions(check_kind="ready")
    runtime_gateway_check = _runtime_gateway_contract_check()
    is_ready = (
        aiida_check.status in {"ok", "skipped"}
        and deep_readiness_check.status in {"ok", "skipped"}
        and real_adapter_check.status in {"ok", "skipped"}
        and runtime_gateway_check["status"] == "ok"
    )
    if not is_ready:
        response.status_code = status.HTTP_503_SERVICE_UNAVAILABLE
    return {
        "status": "ready" if is_ready else "not_ready",
        "checks": {
            "managed_aiida_runtime": aiida_check.as_dict(),
            "user_managed_deep_readiness": deep_readiness_check.as_dict(),
            "slurm_real_adapter_preconditions": real_adapter_check.as_dict(),
            "runtime_gateway_contract": runtime_gateway_check,
        },
    }
