from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request

from app.deps import require_tenant_id, require_user_identity
from app.middleware import get_request_id
from app.schemas.runtime import (
    ExecutionEvent,
    RuntimeEventAck,
    RuntimeJobStatusResponse,
    SubmitJobAccepted,
    SubmitJobCommand,
)
from app.schemas.zpe import ZPEJobStatus
from services.runtime import (
    RuntimeConfigurationError,
    RuntimeConflictError,
    RuntimeDownstreamError,
    RuntimeNotFoundError,
    get_runtime_store,
)

router = APIRouter(prefix="/api/runtime", tags=["runtime"])


@router.post("/jobs:submit", response_model=SubmitJobAccepted, status_code=202)
async def submit_runtime_job(command: SubmitJobCommand, request: Request) -> SubmitJobAccepted:
    user = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    if command.tenant_id != tenant_id:
        raise HTTPException(status_code=403, detail="tenant scope violation")
    if command.requested_by.user_id != user.user_id:
        raise HTTPException(status_code=403, detail="requested_by scope violation")
    trace_id = get_request_id(request)
    try:
        result = get_runtime_store().submit(command, trace_id=trace_id)
    except RuntimeConflictError as exc:
        detail = str(exc)
        status = 409 if detail == "idempotency_conflict" else 422
        raise HTTPException(status_code=status, detail=detail) from exc
    except RuntimeConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except RuntimeDownstreamError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    return result.response


@router.post("/jobs/{job_id}/events", response_model=RuntimeEventAck)
async def post_runtime_event(job_id: str, event: ExecutionEvent, request: Request) -> RuntimeEventAck:
    _ = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    if event.job_id != job_id:
        raise HTTPException(status_code=400, detail="job_id mismatch")
    if event.tenant_id != tenant_id:
        raise HTTPException(status_code=403, detail="tenant scope violation")
    try:
        ack = get_runtime_store().apply_event(event)
    except RuntimeNotFoundError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except RuntimeConflictError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except RuntimeConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except RuntimeDownstreamError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
    return RuntimeEventAck(ok=True, idempotent=ack.idempotent)


@router.get("/jobs/{job_id}", response_model=RuntimeJobStatusResponse)
async def get_runtime_job(job_id: str, request: Request) -> RuntimeJobStatusResponse:
    _ = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    try:
        return get_runtime_store().get_status(job_id, tenant_id=tenant_id)
    except RuntimeNotFoundError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except RuntimeConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except RuntimeDownstreamError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@router.get("/jobs/{job_id}/detail")
async def get_runtime_job_detail(job_id: str, request: Request) -> dict[str, object]:
    _ = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    try:
        return get_runtime_store().get_detail(job_id, tenant_id=tenant_id)
    except RuntimeNotFoundError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except RuntimeConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except RuntimeDownstreamError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc


@router.get("/jobs/{job_id}/projection", response_model=ZPEJobStatus)
async def get_runtime_job_projection(job_id: str, request: Request) -> ZPEJobStatus:
    _ = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    try:
        return get_runtime_store().get_projection_status(job_id, tenant_id=tenant_id)
    except RuntimeNotFoundError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except RuntimeConfigurationError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    except RuntimeDownstreamError as exc:
        raise HTTPException(status_code=502, detail=str(exc)) from exc
