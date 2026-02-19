from __future__ import annotations

from pathlib import Path
from uuid import uuid4

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import PlainTextResponse

from app.deps import require_tenant_id, require_user_identity
from app.middleware import get_request_id
from app.schemas.runtime import (
    ExecutionEvent,
    RuntimeNodeJoinTokenRequest,
    RuntimeNodeJoinTokenResponse,
    RuntimeNodeRegisterRequest,
    RuntimeNodeRegisterResponse,
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
from services.runtime_nodes import get_runtime_node_store
from services.structures import get_structure
from services.zpe.parse import (
    extract_fixed_indices,
    parse_atomic_species,
    parse_kpoints_automatic,
    parse_qe_atoms,
    parse_qe_structure,
)
from services.zpe.settings import get_zpe_settings

router = APIRouter(prefix="/api/runtime", tags=["runtime"])

_REPO_ROOT = Path(__file__).resolve().parents[4]
_COMPUTE_INSTALL_SCRIPT = _REPO_ROOT / "scripts" / "install-compute-node.sh"


def _default_queue_name() -> str:
    configured = get_zpe_settings().queue_name.strip()
    return configured or "zpe"


def _normalize_queue_name(value: str | None) -> str:
    if value is None:
        return _default_queue_name()
    normalized = value.strip()
    return normalized or _default_queue_name()


def _resolve_api_base(request: Request) -> str:
    return str(request.base_url).rstrip("/")


def _build_install_command(
    *,
    install_script_url: str,
    api_base: str,
    join_token: str,
    queue_name: str,
    node_name_hint: str | None,
) -> str:
    command = (
        f"curl -fsSL {install_script_url} | bash -s -- "
        f"--api-base {api_base} --join-token {join_token} --queue-name {queue_name}"
    )
    if node_name_hint and node_name_hint.strip():
        command = f"{command} --name {node_name_hint.strip()}"
    return command


@router.get("/nodes/install.sh", response_class=PlainTextResponse)
async def runtime_compute_node_install_script() -> str:
    if not _COMPUTE_INSTALL_SCRIPT.exists():
        raise HTTPException(status_code=404, detail="install script not found")
    return _COMPUTE_INSTALL_SCRIPT.read_text(encoding="utf-8")


@router.post("/nodes/join-token", response_model=RuntimeNodeJoinTokenResponse)
async def issue_runtime_node_join_token(
    payload: RuntimeNodeJoinTokenRequest,
    request: Request,
) -> RuntimeNodeJoinTokenResponse:
    user = require_user_identity(request)
    tenant_id = require_tenant_id(request)
    queue_name = _normalize_queue_name(payload.queue_name)
    ttl_seconds = payload.ttl_seconds or get_zpe_settings().enroll_token_ttl_seconds
    try:
        token = get_runtime_node_store().create_join_token(
            tenant_id=tenant_id,
            owner_user_id=user.user_id,
            queue_name=queue_name,
            ttl_seconds=ttl_seconds,
            node_name_hint=payload.node_name_hint,
        )
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    api_base = _resolve_api_base(request)
    install_script_url = f"{api_base}/api/runtime/nodes/install.sh"
    return RuntimeNodeJoinTokenResponse(
        join_token=token.token,
        expires_at=token.expires_at,
        token_ttl_seconds=token.ttl_seconds,
        queue_name=queue_name,
        install_script_url=install_script_url,
        register_endpoint=f"{api_base}/api/runtime/nodes/register",
        install_command=_build_install_command(
            install_script_url=install_script_url,
            api_base=api_base,
            join_token=token.token,
            queue_name=queue_name,
            node_name_hint=payload.node_name_hint,
        ),
    )


@router.post("/nodes/register", response_model=RuntimeNodeRegisterResponse)
async def register_runtime_compute_node(
    payload: RuntimeNodeRegisterRequest,
) -> RuntimeNodeRegisterResponse:
    try:
        join_token = get_runtime_node_store().consume_join_token(payload.token)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="join token not found") from exc
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc

    server_id = f"compute-{uuid4().hex}"
    queue_name = _normalize_queue_name(join_token.queue_name)
    target_store = get_runtime_node_store()
    target = target_store.add_target(
        tenant_id=join_token.tenant_id,
        user_id=join_token.owner_user_id,
        queue_name=queue_name,
        server_id=server_id,
        name=payload.name,
        metadata=payload.meta,
    )
    if target_store.get_active_target(join_token.tenant_id, join_token.owner_user_id) is None:
        target_store.set_active_target(
            join_token.tenant_id,
            join_token.owner_user_id,
            target.target_id,
        )

    return RuntimeNodeRegisterResponse(
        server_id=server_id,
        target_id=target.target_id,
        queue_name=target.queue_name,
        registered_at=target.registered_at,
        name=target.name,
    )


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
    tenant_id = require_tenant_id(raw)
    try:
        target_store = get_runtime_node_store()
        targets = target_store.list_targets(tenant_id, user.user_id)
        active = target_store.get_active_target(tenant_id, user.user_id)
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    total = len(targets)
    sliced = targets[offset : offset + limit]
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
    tenant_id = require_tenant_id(raw)
    try:
        target_store = get_runtime_node_store()
        target_store.ensure_target_owner(tenant_id, user.user_id, target_id)
        target_store.set_active_target(tenant_id, user.user_id, target_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="target not found") from exc
    except RuntimeError as exc:
        raise HTTPException(status_code=503, detail=str(exc)) from exc
    return QueueTargetSelectResponse(active_target_id=target_id)
