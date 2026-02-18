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
from app.schemas.zpe import (
    QueueTarget,
    QueueTargetListResponse,
    QueueTargetSelectResponse,
    ZPEJobStatus,
    ZPEParseRequest,
    ZPEParseResponse,
)
from app.schemas.common import Pagination
from services.runtime import (
    RuntimeConfigurationError,
    RuntimeConflictError,
    RuntimeDownstreamError,
    RuntimeNotFoundError,
    get_runtime_store,
)
from services.structures import get_structure
from services.zpe.parse import (
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_qe_atoms,
    parse_qe_structure,
)
from services.zpe.queue_targets import get_queue_target_store

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


@router.post("/parse", response_model=ZPEParseResponse)
async def parse_runtime_input(request: ZPEParseRequest) -> ZPEParseResponse:
    if request.structure_id:
        try:
            structure = get_structure(request.structure_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail="Structure not found") from exc
        fixed_indices = extract_fixed_indices(request.content)
        atoms = parse_qe_atoms(request.content)
        if len(atoms) != len(structure.atoms):
            raise HTTPException(
                status_code=409,
                detail="Structure does not match QE input atom count",
            )
    else:
        structure, fixed_indices = parse_qe_structure(request.content)
    kpts = parse_kpoints_automatic(request.content)
    return ZPEParseResponse(
        structure=structure,
        fixed_indices=fixed_indices,
        atomic_species=parse_atomic_species(request.content),
        kpoints=kpts[0] if kpts else None,
    )


@router.get("/targets", response_model=QueueTargetListResponse)
async def list_runtime_targets(
    raw: Request,
    limit: int = 50,
    offset: int = 0,
) -> QueueTargetListResponse:
    user = require_user_identity(raw)
    target_store = get_queue_target_store()
    targets = target_store.list_targets(user.user_id)
    total = len(targets)
    sliced = targets[offset : offset + limit]
    active = target_store.get_active_target(user.user_id)
    response_targets = [
        QueueTarget(
            id=target.target_id,
            queue_name=target.queue_name,
            server_id=target.server_id,
            registered_at=target.registered_at,
            name=target.name,
        )
        for target in sliced
    ]
    return QueueTargetListResponse(
        targets=response_targets,
        active_target_id=active.target_id if active else None,
        pagination=Pagination(total=total, limit=limit, offset=offset),
    )


@router.put("/targets/{target_id}/active", response_model=QueueTargetSelectResponse)
async def select_runtime_target(
    target_id: str,
    raw: Request,
) -> QueueTargetSelectResponse:
    user = require_user_identity(raw)
    target_store = get_queue_target_store()
    try:
        target_store.ensure_target_owner(user.user_id, target_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="target not found") from exc
    target_store.set_active_target(user.user_id, target_id)
    return QueueTargetSelectResponse(active_target_id=target_id)
