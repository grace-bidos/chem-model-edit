from __future__ import annotations

import logging

from fastapi import APIRouter, HTTPException, Query, Request

from app.deps import require_admin, require_user_identity
from app.middleware import get_request_id
from app.schemas.common import Pagination
from app.schemas.zpe import (
    OpsFlagsRequest,
    OpsFlagsResponse,
    QueueTarget,
    QueueTargetListResponse,
    QueueTargetSelectResponse,
    ZPEParseRequest,
    ZPEParseResponse,
)
from services.structures import get_structure
from services.zpe import get_zpe_settings
from services.zpe.ops_flags import get_ops_flags, set_ops_flags
from services.zpe.parse import (
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_qe_atoms,
    parse_qe_structure,
)
from services.zpe.queue_targets import get_queue_target_store
from services.zpe.structured_log import log_event

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/zpe", tags=["zpe"])


@router.post("/parse", response_model=ZPEParseResponse)
async def zpe_parse(request: ZPEParseRequest) -> ZPEParseResponse:
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
async def list_queue_targets(
    raw: Request,
    limit: int = Query(50, ge=1),
    offset: int = Query(0, ge=0),
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
async def select_queue_target(
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


@router.get("/admin/ops", response_model=OpsFlagsResponse)
async def zpe_ops_flags(raw: Request) -> OpsFlagsResponse:
    require_admin(raw)
    flags = get_ops_flags()
    return OpsFlagsResponse(
        submission_enabled=flags.submission_enabled,
        dequeue_enabled=flags.dequeue_enabled,
        result_read_source=flags.result_read_source,
    )


@router.patch("/admin/ops", response_model=OpsFlagsResponse)
async def zpe_ops_flags_update(
    request: OpsFlagsRequest,
    raw: Request,
) -> OpsFlagsResponse:
    require_admin(raw)
    flags = set_ops_flags(
        submission_enabled=request.submission_enabled,
        dequeue_enabled=request.dequeue_enabled,
        result_read_source=request.result_read_source,
    )
    settings = get_zpe_settings()
    log_event(
        logger,
        event="zpe_ops_flags_updated",
        service="control-plane",
        stage="ops",
        status="updated",
        request_id=get_request_id(raw),
        submission_enabled=flags.submission_enabled,
        dequeue_enabled=flags.dequeue_enabled,
        result_read_source=flags.result_read_source,
        backend=settings.compute_mode,
        result_store=settings.result_store,
    )
    return OpsFlagsResponse(
        submission_enabled=flags.submission_enabled,
        dequeue_enabled=flags.dequeue_enabled,
        result_read_source=flags.result_read_source,
    )
