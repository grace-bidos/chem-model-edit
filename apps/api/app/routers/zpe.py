from __future__ import annotations

import logging
from typing import Any, Literal

from fastapi import APIRouter, HTTPException, Query, Request
from fastapi.responses import Response
from pydantic import ValidationError

from app.deps import (
    has_admin_access,
    require_admin,
    require_job_owner,
    require_tenant_id,
    require_user_identity,
    require_worker,
)
from app.middleware import get_request_id, get_tenant_id
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
from services.zpe import get_zpe_settings
from services.zpe.backends import (
    IdempotencyKeyConflictError,
    JobConflictError,
    NoActiveQueueTargetError,
    ResultArtifactMissingError,
    StructureMismatchError,
    get_job_file,
    get_job_result,
    get_job_status,
    submit_compute_failure,
    submit_compute_result,
    submit_job,
)
from services.zpe.enroll import get_enroll_store
from services.zpe.job_meta import get_job_meta_store
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
from services.zpe.slurm_policy import (
    SlurmPolicyConfigError,
    SlurmPolicyDeniedError,
)
from services.zpe.structured_log import log_event, write_audit_event
from services.zpe.worker_auth import get_worker_token_store

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/zpe", tags=["zpe"])


def _require_legacy_worker_endpoints(flags: OpsFlags) -> None:
    if not flags.legacy_worker_endpoints_enabled:
        raise HTTPException(
            status_code=503,
            detail="legacy worker endpoints disabled by cutover flag",
        )


def _require_job_tenant_meta(job_id: str) -> dict[str, Any]:
    meta_store = get_job_meta_store()
    job_meta = meta_store.get_meta(job_id)
    tenant_id = job_meta.get("tenant_id")
    if not isinstance(tenant_id, str) or not tenant_id.strip():
        raise HTTPException(status_code=409, detail="job tenant metadata missing")
    return job_meta


def _write_audit_or_raise_unavailable(
    *,
    trace_request_id: str,
    event: str,
    **fields: Any,
) -> None:
    try:
        write_audit_event(**fields)
    except Exception as exc:
        log_event(
            logger,
            event="zpe_audit_write_failed",
            service="control-plane",
            stage="audit",
            status="error",
            request_id=trace_request_id,
            source_event=event,
            error_message=str(exc),
        )
        raise HTTPException(status_code=503, detail="audit sink unavailable") from exc


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
    tenant_id = require_tenant_id(raw)
    request_id = get_request_id(raw)
    settings = get_zpe_settings()
    flags = get_ops_flags()
    if not flags.submission_enabled:
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
    try:
        outcome = submit_job(
            request,
            user_id=user.user_id,
            tenant_id=tenant_id,
            request_id=request_id,
        )
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Structure not found") from exc
    except NoActiveQueueTargetError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except StructureMismatchError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except SlurmPolicyDeniedError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except IdempotencyKeyConflictError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except SlurmPolicyConfigError as exc:
        raise HTTPException(
            status_code=503,
            detail=f"slurm policy configuration error: {exc}",
        ) from exc
    resolution = outcome.slurm_resolution
    log_event(
        logger,
        event="zpe_job_enqueued",
        service="control-plane",
        stage="enqueue",
        status="queued",
        request_id=request_id,
        tenant_id=tenant_id,
        job_id=outcome.job_id,
        user_id=user.user_id,
        structure_id=request.structure_id,
        queue_name=outcome.resolved_queue_name,
        requested_queue_name=outcome.requested_queue_name,
        queue_resolution_used_fallback=(
            resolution.used_fallback if resolution is not None else False
        ),
        slurm_partition=resolution.mapping.partition if resolution else None,
        slurm_account=resolution.mapping.account if resolution else None,
        slurm_qos=resolution.mapping.qos if resolution else None,
        backend=settings.compute_mode,
        result_store=settings.result_store,
    )
    _write_audit_or_raise_unavailable(
        trace_request_id=request_id,
        event="zpe_job_enqueued",
        event_id=f"zpe_job_enqueued:{outcome.job_id}:{request_id}",
        request_id=request_id,
        tenant_id=tenant_id,
        actor_id=user.user_id,
        operation="execute",
        resource_type="zpe_job",
        resource_id=outcome.job_id,
        outcome="succeeded",
        metadata={
            "calc_mode": request.calc_mode,
            "queue_name": outcome.resolved_queue_name,
        },
    )
    return ZPEJobResponse(id=outcome.job_id)


@router.get("/jobs/{job_id}", response_model=ZPEJobStatus)
async def zpe_job_status(job_id: str, raw: Request) -> ZPEJobStatus:
    require_job_owner(raw, job_id)
    try:
        status = get_job_status(job_id)
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
    try:
        result_dict = get_job_result(job_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except JobConflictError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ResultArtifactMissingError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
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
    try:
        payload = get_job_file(job_id, kind)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="job not found") from exc
    except JobConflictError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except ResultArtifactMissingError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
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
    lease_tenant_id = meta.get("tenant_id")
    if isinstance(lease_tenant_id, str) and lease_tenant_id.strip():
        raw.state.tenant_id = lease_tenant_id.strip()
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
        tenant_id=get_tenant_id(raw),
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
    job_meta = _require_job_tenant_meta(job_id)
    tenant_id = str(job_meta["tenant_id"])
    if request.tenant_id != tenant_id:
        raise HTTPException(status_code=403, detail="tenant boundary violation")
    execution_event = request.execution_event
    if execution_event is not None:
        if execution_event.tenant_id != tenant_id:
            raise HTTPException(status_code=403, detail="tenant boundary violation")
        if execution_event.job_id != job_id:
            raise HTTPException(status_code=409, detail="execution_event job_id mismatch")
    raw.state.tenant_id = tenant_id
    try:
        submission = submit_compute_result(
            job_id=job_id,
            worker_id=worker_id,
            lease_id=request.lease_id,
            event_id=(execution_event.event_id if execution_event else None),
            result=request.result,
            summary_text=request.summary_text,
            freqs_csv=request.freqs_csv,
            worker_meta=request.meta,
        )
    except PermissionError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    settings = get_zpe_settings()
    job_meta = submission.job_meta
    trace_id = (
        execution_event.trace_id
        if execution_event is not None
        else str(job_meta.get("request_id", get_request_id(raw)))
    )
    log_event(
        logger,
        event="zpe_result_submitted",
        service="control-plane",
        stage="result",
        status="submitted",
        request_id=trace_id,
        trace_id=trace_id,
        tenant_id=tenant_id,
        workspace_id=(
            execution_event.workspace_id
            if execution_event is not None
            else str(job_meta.get("workspace_id", tenant_id))
        ),
        submission_id=(
            execution_event.submission_id
            if execution_event is not None
            else str(job_meta.get("submission_id", job_meta.get("request_id", "")))
        ),
        execution_id=(
            execution_event.execution_id
            if execution_event is not None
            else str(job_meta.get("execution_id", request.lease_id))
        ),
        event_id=(
            execution_event.event_id
            if execution_event is not None
            else f"zpe_result_submitted:{job_id}:{request.lease_id}"
        ),
        state=(execution_event.state if execution_event is not None else "completed"),
        job_id=job_id,
        worker_id=worker_id,
        user_id=job_meta.get("user_id"),
        slurm_job_id=(
            execution_event.scheduler_ref.slurm_job_id
            if execution_event is not None and execution_event.scheduler_ref is not None
            else None
        ),
        slurm_partition=(
            execution_event.scheduler_ref.partition
            if execution_event is not None and execution_event.scheduler_ref is not None
            else None
        ),
        slurm_qos=(
            execution_event.scheduler_ref.qos
            if execution_event is not None and execution_event.scheduler_ref is not None
            else None
        ),
        backend=settings.compute_mode,
        result_store=settings.result_store,
        duration_ms=submission.duration_ms,
        exit_code=0,
        qe_version=submission.qe_version,
    )
    outcome = submission.outcome
    request_id = trace_id
    _write_audit_or_raise_unavailable(
        trace_request_id=request_id,
        event="zpe_result_submitted",
        event_id=(
            execution_event.event_id
            if execution_event is not None
            else f"zpe_result_submitted:{job_id}:{request.lease_id}"
        ),
        request_id=request_id,
        tenant_id=tenant_id,
        actor_id=worker_id,
        operation="update",
        resource_type="zpe_job",
        resource_id=job_id,
        outcome="succeeded",
        metadata={"idempotent": outcome.idempotent, "endpoint": "result"},
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
    job_meta = _require_job_tenant_meta(job_id)
    tenant_id = str(job_meta["tenant_id"])
    if request.tenant_id != tenant_id:
        raise HTTPException(status_code=403, detail="tenant boundary violation")
    execution_event = request.execution_event
    if execution_event is not None:
        if execution_event.tenant_id != tenant_id:
            raise HTTPException(status_code=403, detail="tenant boundary violation")
        if execution_event.job_id != job_id:
            raise HTTPException(status_code=409, detail="execution_event job_id mismatch")
        if execution_event.error.code != request.error_code:
            raise HTTPException(
                status_code=409,
                detail="execution_event error.code mismatch",
            )
        if execution_event.error.message != request.error_message:
            raise HTTPException(
                status_code=409,
                detail="execution_event error.message mismatch",
            )
    raw.state.tenant_id = tenant_id
    try:
        submission = submit_compute_failure(
            job_id=job_id,
            worker_id=worker_id,
            lease_id=request.lease_id,
            event_id=(execution_event.event_id if execution_event else None),
            error_code=request.error_code,
            error_message=request.error_message,
            traceback=request.traceback,
        )
    except PermissionError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    settings = get_zpe_settings()
    job_meta = submission.job_meta
    outcome = submission.outcome
    trace_id = (
        execution_event.trace_id
        if execution_event is not None
        else str(job_meta.get("request_id", get_request_id(raw)))
    )
    log_event(
        logger,
        event="zpe_result_failed",
        service="control-plane",
        stage="result",
        status="failed",
        request_id=trace_id,
        trace_id=trace_id,
        tenant_id=tenant_id,
        workspace_id=(
            execution_event.workspace_id
            if execution_event is not None
            else str(job_meta.get("workspace_id", tenant_id))
        ),
        submission_id=(
            execution_event.submission_id
            if execution_event is not None
            else str(job_meta.get("submission_id", job_meta.get("request_id", "")))
        ),
        execution_id=(
            execution_event.execution_id
            if execution_event is not None
            else str(job_meta.get("execution_id", request.lease_id))
        ),
        event_id=(
            execution_event.event_id
            if execution_event is not None
            else f"zpe_result_failed:{job_id}:{request.lease_id}"
        ),
        state=(execution_event.state if execution_event is not None else "failed"),
        job_id=job_id,
        worker_id=worker_id,
        user_id=job_meta.get("user_id"),
        backend=settings.compute_mode,
        result_store=settings.result_store,
        exit_code=1,
        qe_version=None,
        error_code=request.error_code,
        error_retryable=(
            execution_event.error.retryable if execution_event is not None else None
        ),
        retry_count=outcome.retry_count,
        requeued=outcome.requeued,
    )
    request_id = trace_id
    _write_audit_or_raise_unavailable(
        trace_request_id=request_id,
        event="zpe_result_failed",
        event_id=(
            execution_event.event_id
            if execution_event is not None
            else f"zpe_result_failed:{job_id}:{request.lease_id}"
        ),
        request_id=request_id,
        tenant_id=tenant_id,
        actor_id=worker_id,
        operation="update",
        resource_type="zpe_job",
        resource_id=job_id,
        outcome="succeeded",
        metadata={
            "error_code": request.error_code,
            "retry_count": outcome.retry_count,
            "requeued": outcome.requeued,
            "endpoint": "failed",
        },
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
