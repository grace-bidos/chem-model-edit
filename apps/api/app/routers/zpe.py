from __future__ import annotations

import logging
from typing import Literal

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import Response
from pydantic import ValidationError

from app.deps import (
    has_admin_access,
    require_admin,
    require_job_owner,
    require_user_identity,
    require_worker,
)
from app.middleware import get_request_id
from app.schemas.common import Pagination
from app.schemas.zpe import (
    ComputeFailedRequest,
    ComputeFailedResponse,
    ComputeLeaseResponse,
    ComputeRegisterRequest,
    ComputeRegisterResponse,
    ComputeResultRequest,
    ComputeResultResponse,
    ComputeRevokeResponse,
    EnrollTokenRequest,
    EnrollTokenResponse,
    OpsFlagsRequest,
    OpsFlagsResponse,
    QueueTarget,
    QueueTargetListResponse,
    QueueTargetSelectResponse,
    ZPEJobRequest,
    ZPEJobResponse,
    ZPEJobResultResponse,
    ZPEJobStatus,
    ZPEParseRequest,
    ZPEParseResponse,
    ZPEResult,
)
from services.structures import get_structure
from services.zpe import (
    ensure_mobile_indices,
    enqueue_zpe_job,
    get_result_store,
    get_zpe_settings,
)
from services.zpe.compute_results import submit_failure, submit_result
from services.zpe.enroll import get_enroll_store
from services.zpe.job_meta import get_job_meta_store
from services.zpe.job_owner import get_job_owner_store
from services.zpe.lease import lease_next_job
from services.zpe.ops_flags import OpsFlags, get_ops_flags, set_ops_flags
from services.zpe.parse import (
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_qe_atoms,
    parse_qe_structure,
)
from services.zpe.queue_targets import get_queue_target_store
from services.zpe.structured_log import log_event
from services.zpe.worker_auth import get_worker_token_store

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/zpe", tags=["zpe"])


def _ensure_job_finished(job_id: str) -> None:
    store = get_result_store()
    try:
        status = store.get_status(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    if status.status == "failed":
        detail = status.detail or "unknown"
        raise HTTPException(status_code=409, detail=f"job failed: {detail}")
    if status.status != "finished":
        raise HTTPException(
            status_code=409,
            detail=f"job not finished (status={status.status})",
        )


def _require_legacy_submission_route(flags: OpsFlags) -> None:
    if flags.submission_route != "redis-worker":
        raise HTTPException(status_code=503, detail="submission route not available")


def _require_legacy_result_read_source(flags: OpsFlags) -> None:
    if flags.result_read_source != "redis":
        raise HTTPException(status_code=503, detail="result read source not available")


def _require_legacy_worker_endpoints(flags: OpsFlags) -> None:
    if not flags.legacy_worker_endpoints_enabled:
        raise HTTPException(
            status_code=503,
            detail="legacy worker endpoints disabled by cutover flag",
        )


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


@router.post("/jobs", response_model=ZPEJobResponse)
async def zpe_jobs(request: ZPEJobRequest, raw: Request) -> ZPEJobResponse:
    user = require_user_identity(raw)
    request_id = get_request_id(raw)
    flags = get_ops_flags()
    _require_legacy_submission_route(flags)
    if not flags.submission_enabled:
        settings = get_zpe_settings()
        log_event(
            logger,
            event="zpe_submission_blocked",
            service="control-plane",
            stage="enqueue",
            status="blocked",
            request_id=request_id,
            user_id=user.user_id,
            backend=settings.compute_mode,
            result_store=settings.result_store,
        )
        raise HTTPException(status_code=503, detail="zpe submissions disabled")
    target_store = get_queue_target_store()
    target = target_store.get_active_target(user.user_id)
    if not target:
        raise HTTPException(status_code=400, detail="no active queue target")
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
    mobile_indices = ensure_mobile_indices(
        request.mobile_indices, len(structure.atoms), fixed_indices
    )
    payload = request.model_dump()
    payload["mobile_indices"] = mobile_indices
    job_id = enqueue_zpe_job(payload, queue_name=target.queue_name)
    owner_store = get_job_owner_store()
    owner_store.set_owner(job_id, user.user_id)
    meta_store = get_job_meta_store()
    meta_store.set_meta(
        job_id,
        {
            "request_id": request_id,
            "user_id": user.user_id,
            "structure_id": request.structure_id,
            "queue_name": target.queue_name,
            "calc_mode": request.calc_mode,
        },
    )
    settings = get_zpe_settings()
    log_event(
        logger,
        event="zpe_job_enqueued",
        service="control-plane",
        stage="enqueue",
        status="queued",
        request_id=request_id,
        job_id=job_id,
        user_id=user.user_id,
        structure_id=request.structure_id,
        queue_name=target.queue_name,
        backend=settings.compute_mode,
        result_store=settings.result_store,
    )
    return ZPEJobResponse(id=job_id)


@router.get("/jobs/{job_id}", response_model=ZPEJobStatus)
async def zpe_job_status(job_id: str, raw: Request) -> ZPEJobStatus:
    require_job_owner(raw, job_id)
    flags = get_ops_flags()
    _require_legacy_result_read_source(flags)
    store = get_result_store()
    try:
        status = store.get_status(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    return ZPEJobStatus(
        status=status.status,
        detail=status.detail,
        updated_at=status.updated_at,
    )


@router.get("/jobs/{job_id}/result", response_model=ZPEJobResultResponse)
async def zpe_job_result(job_id: str, raw: Request) -> ZPEJobResultResponse:
    require_job_owner(raw, job_id)
    flags = get_ops_flags()
    _require_legacy_result_read_source(flags)
    store = get_result_store()
    _ensure_job_finished(job_id)
    try:
        result_dict = store.get_result(job_id)
    except KeyError as exc:
        raise HTTPException(
            status_code=500, detail="result missing after completion"
        ) from exc
    if "kpts" in result_dict and "kpoints" not in result_dict:
        result_dict = {
            **{key: value for key, value in result_dict.items() if key != "kpts"},
            "kpoints": result_dict.get("kpts"),
        }
    try:
        result = ZPEResult(**result_dict)
    except ValidationError as exc:
        raise HTTPException(status_code=500, detail="result data invalid") from exc
    return ZPEJobResultResponse(result=result)


@router.get("/jobs/{job_id}/files")
async def zpe_job_files(
    job_id: str,
    raw: Request,
    kind: Literal["summary", "freqs"] = Query(...),
) -> Response:
    require_job_owner(raw, job_id)
    flags = get_ops_flags()
    _require_legacy_result_read_source(flags)
    store = get_result_store()
    _ensure_job_finished(job_id)
    try:
        payload = store.get_file(job_id, kind)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except KeyError as exc:
        raise HTTPException(
            status_code=500, detail="file missing after completion"
        ) from exc
    filename = "summary.txt" if kind == "summary" else "freqs.csv"
    media_type = "text/plain" if kind == "summary" else "text/csv"
    return Response(
        content=payload,
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@router.post("/compute/enroll-tokens", response_model=EnrollTokenResponse)
async def zpe_compute_enroll_token(
    request: EnrollTokenRequest,
    raw: Request,
) -> EnrollTokenResponse:
    flags = get_ops_flags()
    _require_legacy_worker_endpoints(flags)
    owner_id = None
    if not has_admin_access(raw):
        user = require_user_identity(raw)
        owner_id = user.user_id
    if request.ttl_seconds is not None and request.ttl_seconds <= 0:
        raise HTTPException(status_code=400, detail="ttl_seconds must be >= 1")
    store = get_enroll_store()
    token = store.create_token(
        ttl_seconds=request.ttl_seconds,
        label=request.label,
        owner_id=owner_id,
    )
    return EnrollTokenResponse(
        token=token.token,
        expires_at=token.expires_at,
        ttl_seconds=token.ttl_seconds,
        label=token.label,
    )


@router.post("/compute/servers", response_model=ComputeRegisterResponse)
async def zpe_compute_register(
    request: ComputeRegisterRequest,
) -> ComputeRegisterResponse:
    flags = get_ops_flags()
    _require_legacy_worker_endpoints(flags)
    store = get_enroll_store()
    try:
        registration = store.consume_token(
            request.token,
            name=request.name,
            meta=request.meta,
        )
    except KeyError as exc:
        raise HTTPException(status_code=400, detail="invalid enroll token") from exc
    queue_name = request.queue_name or get_zpe_settings().queue_name
    if registration.owner_id:
        target_store = get_queue_target_store()
        had_active_target = (
            target_store.get_active_target(registration.owner_id) is not None
        )
        target = target_store.add_target(
            user_id=registration.owner_id,
            queue_name=queue_name,
            server_id=registration.server_id,
            name=request.name,
        )
        if request.activate_target or not had_active_target:
            target_store.set_active_target(registration.owner_id, target.target_id)
    token_store = get_worker_token_store()
    worker_token = token_store.create_token(registration.server_id, label=request.name)
    return ComputeRegisterResponse(
        id=registration.server_id,
        registered_at=registration.registered_at,
        name=registration.name,
        worker_token=worker_token.token,
        token_expires_at=worker_token.expires_at,
        token_ttl_seconds=worker_token.ttl_seconds,
    )


@router.delete("/compute/servers/{server_id}", response_model=ComputeRevokeResponse)
async def zpe_compute_revoke(
    server_id: str,
    raw: Request,
) -> ComputeRevokeResponse:
    require_admin(raw)
    token_store = get_worker_token_store()
    revoked = token_store.revoke_tokens_for_worker(server_id)
    return ComputeRevokeResponse(revoked_count=revoked)


@router.post("/compute/jobs/lease", response_model=ComputeLeaseResponse)
async def zpe_compute_lease(raw: Request) -> Response | ComputeLeaseResponse:
    worker_id = require_worker(raw)
    request_id = get_request_id(raw)
    flags = get_ops_flags()
    _require_legacy_worker_endpoints(flags)
    if not flags.dequeue_enabled:
        settings = get_zpe_settings()
        log_event(
            logger,
            event="zpe_dequeue_paused",
            service="control-plane",
            stage="lease",
            status="paused",
            request_id=request_id,
            worker_id=worker_id,
            backend=settings.compute_mode,
            result_store=settings.result_store,
        )
        return Response(status_code=204)
    lease = lease_next_job(worker_id)
    if lease is None:
        return Response(status_code=204)
    meta = lease.meta
    try:
        payload = ZPEJobRequest(**lease.payload).model_dump(
            by_alias=True, exclude_none=True
        )
    except ValidationError:
        payload = lease.payload
    settings = get_zpe_settings()
    log_event(
        logger,
        event="zpe_lease_granted",
        service="control-plane",
        stage="lease",
        status="granted",
        request_id=meta.get("request_id", request_id),
        job_id=lease.job_id,
        worker_id=worker_id,
        user_id=meta.get("user_id"),
        backend=settings.compute_mode,
        result_store=settings.result_store,
    )
    return ComputeLeaseResponse(
        job_id=lease.job_id,
        payload=payload,
        lease_id=lease.lease_id,
        lease_ttl_seconds=lease.lease_ttl_seconds,
        meta=meta,
    )


@router.post("/compute/jobs/{job_id}/result", response_model=ComputeResultResponse)
async def zpe_compute_result(
    job_id: str,
    request: ComputeResultRequest,
    raw: Request,
) -> ComputeResultResponse:
    worker_id = require_worker(raw)
    flags = get_ops_flags()
    _require_legacy_worker_endpoints(flags)
    meta_store = get_job_meta_store()
    job_meta = meta_store.get_meta(job_id)
    try:
        outcome = submit_result(
            job_id=job_id,
            worker_id=worker_id,
            lease_id=request.lease_id,
            result=request.result,
            summary_text=request.summary_text,
            freqs_csv=request.freqs_csv,
        )
    except PermissionError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    settings = get_zpe_settings()
    duration_ms = None
    qe_version = request.meta.get("qe_version") if request.meta else None
    if request.meta and "computation_time_seconds" in request.meta:
        try:
            duration_ms = int(float(request.meta["computation_time_seconds"]) * 1000)
        except (TypeError, ValueError):
            duration_ms = None
    log_event(
        logger,
        event="zpe_result_submitted",
        service="control-plane",
        stage="result",
        status="submitted",
        request_id=job_meta.get("request_id", get_request_id(raw)),
        job_id=job_id,
        worker_id=worker_id,
        user_id=job_meta.get("user_id"),
        backend=settings.compute_mode,
        result_store=settings.result_store,
        duration_ms=duration_ms,
        exit_code=0,
        qe_version=qe_version,
    )
    return ComputeResultResponse(idempotent=outcome.idempotent)


@router.post("/compute/jobs/{job_id}/failed", response_model=ComputeFailedResponse)
async def zpe_compute_failed(
    job_id: str,
    request: ComputeFailedRequest,
    raw: Request,
) -> ComputeFailedResponse:
    worker_id = require_worker(raw)
    flags = get_ops_flags()
    _require_legacy_worker_endpoints(flags)
    meta_store = get_job_meta_store()
    job_meta = meta_store.get_meta(job_id)
    try:
        outcome = submit_failure(
            job_id=job_id,
            worker_id=worker_id,
            lease_id=request.lease_id,
            error_code=request.error_code,
            error_message=request.error_message,
            traceback=request.traceback,
        )
    except PermissionError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    settings = get_zpe_settings()
    log_event(
        logger,
        event="zpe_result_failed",
        service="control-plane",
        stage="result",
        status="failed",
        request_id=job_meta.get("request_id", get_request_id(raw)),
        job_id=job_id,
        worker_id=worker_id,
        user_id=job_meta.get("user_id"),
        backend=settings.compute_mode,
        result_store=settings.result_store,
        exit_code=1,
        qe_version=None,
        error_code=request.error_code,
        retry_count=outcome.retry_count,
        requeued=outcome.requeued,
    )
    return ComputeFailedResponse(
        requeued=outcome.requeued,
        retry_count=outcome.retry_count,
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
        submission_route=flags.submission_route,
        result_read_source=flags.result_read_source,
        legacy_worker_endpoints_enabled=flags.legacy_worker_endpoints_enabled,
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
        submission_route=request.submission_route,
        result_read_source=request.result_read_source,
        legacy_worker_endpoints_enabled=request.legacy_worker_endpoints_enabled,
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
        submission_route=flags.submission_route,
        result_read_source=flags.result_read_source,
        legacy_worker_endpoints_enabled=flags.legacy_worker_endpoints_enabled,
        backend=settings.compute_mode,
        result_store=settings.result_store,
    )
    return OpsFlagsResponse(
        submission_enabled=flags.submission_enabled,
        dequeue_enabled=flags.dequeue_enabled,
        submission_route=flags.submission_route,
        result_read_source=flags.result_read_source,
        legacy_worker_endpoints_enabled=flags.legacy_worker_endpoints_enabled,
    )
