from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict
from uuid import uuid4

from app.schemas.zpe import ZPEJobRequest
from services.structures import get_structure
from .io import format_freqs_csv, format_summary
from .job_meta import get_job_meta_store
from .parse import (
    ensure_mobile_indices,
    extract_fixed_indices,
    parse_kpoints_automatic,
    parse_qe_atoms,
    parse_qe_structure,
)
from .http_queue import enqueue_http_job
from .queue import get_queue
from .queue_targets import get_queue_target_store
from .result_store import ResultStore, ZPEJobStatus, get_result_store
from .settings import ZPESettings, get_zpe_settings
from .slurm_policy import (
    SlurmQueueResolution,
    resolve_runtime_slurm_queue,
)
from .submit_idempotency import (
    SubmitIdempotencyRecord,
    compute_submit_request_fingerprint,
    get_submit_idempotency_store,
)
from .thermo import calc_zpe_and_s_vib, normalize_frequencies
from .compute_results import (
    FailureOutcome,
    ResultSubmitOutcome,
    submit_failure,
    submit_result,
)


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _default_kpts() -> tuple[int, int, int]:
    return (1, 1, 1)


class NoActiveQueueTargetError(RuntimeError):
    """Raised when the authenticated user has no active queue target."""


class JobConflictError(RuntimeError):
    """Raised when result/file reads are requested before job completion."""


class ResultArtifactMissingError(RuntimeError):
    """Raised when terminal job artifacts are unexpectedly missing."""


class StructureMismatchError(ValueError):
    """Raised when submitted QE input does not match the requested structure."""


class IdempotencyKeyConflictError(ValueError):
    """Raised when a request id is reused for a different submission payload."""


@dataclass(frozen=True)
class SubmitJobOutcome:
    job_id: str
    requested_queue_name: str
    resolved_queue_name: str
    slurm_resolution: SlurmQueueResolution | None
    idempotent: bool = False


@dataclass(frozen=True)
class SubmitComputeResultOutcome:
    outcome: ResultSubmitOutcome
    job_meta: Dict[str, Any]
    duration_ms: int | None
    qe_version: str | None


@dataclass(frozen=True)
class SubmitComputeFailureOutcome:
    outcome: FailureOutcome
    job_meta: Dict[str, Any]


def _resolve_submission_queue(
    requested_queue: str, *, settings: ZPESettings
) -> SlurmQueueResolution | None:
    if not settings.slurm_policy_path:
        return None
    return resolve_runtime_slurm_queue(
        requested_queue,
        policy_path=Path(settings.slurm_policy_path),
    )


def submit_job(
    request: ZPEJobRequest, *, user_id: str, tenant_id: str, request_id: str
) -> SubmitJobOutcome:
    settings = get_zpe_settings()
    idempotency_store = get_submit_idempotency_store()
    target_store = get_queue_target_store()
    target = target_store.get_active_target(user_id)
    if not target:
        raise NoActiveQueueTargetError("no active queue target")
    resolution = _resolve_submission_queue(target.queue_name, settings=settings)
    resolved_queue_name = resolution.resolved_queue if resolution else target.queue_name
    if request.structure_id:
        structure = get_structure(request.structure_id)
        fixed_indices = extract_fixed_indices(request.content)
        atoms = parse_qe_atoms(request.content)
        if len(atoms) != len(structure.atoms):
            raise StructureMismatchError("Structure does not match QE input atom count")
    else:
        structure, fixed_indices = parse_qe_structure(request.content)
    mobile_indices = ensure_mobile_indices(
        request.mobile_indices, len(structure.atoms), fixed_indices
    )
    payload = request.model_dump()
    payload["mobile_indices"] = mobile_indices
    request_fingerprint = compute_submit_request_fingerprint(
        {
            "request": payload,
            "management_node_id": target.server_id,
            "requested_queue_name": target.queue_name,
            "requested_by_user_id": user_id,
        }
    )
    try:
        claim = idempotency_store.claim_or_get(
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
            request_fingerprint=request_fingerprint,
        )
    except ValueError as exc:
        raise IdempotencyKeyConflictError(str(exc)) from exc
    if claim.state == "ready":
        if claim.record is None:
            raise RuntimeError("submit idempotency record missing")
        return SubmitJobOutcome(
            job_id=claim.record.job_id,
            requested_queue_name=claim.record.requested_queue_name,
            resolved_queue_name=claim.record.resolved_queue_name,
            slurm_resolution=None,
            idempotent=True,
        )
    if claim.state == "pending":
        try:
            existing = idempotency_store.wait_for_record(
                user_id=user_id,
                tenant_id=tenant_id,
                request_id=request_id,
                request_fingerprint=request_fingerprint,
            )
        except ValueError as exc:
            raise IdempotencyKeyConflictError(str(exc)) from exc
        if existing is None:
            raise RuntimeError("idempotency request still in progress; retry")
        return SubmitJobOutcome(
            job_id=existing.job_id,
            requested_queue_name=existing.requested_queue_name,
            resolved_queue_name=existing.resolved_queue_name,
            slurm_resolution=None,
            idempotent=True,
        )
    if claim.claim_token is None:
        raise RuntimeError("submit idempotency claim token missing")
    try:
        job_id = enqueue_zpe_job(payload, queue_name=resolved_queue_name)
    except Exception:
        idempotency_store.release_claim(
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
            claim_token=claim.claim_token,
        )
        raise
    idempotency_store.finalize_claim(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        claim_token=claim.claim_token,
        record=SubmitIdempotencyRecord(
            job_id=job_id,
            request_fingerprint=request_fingerprint,
            requested_queue_name=target.queue_name,
            resolved_queue_name=resolved_queue_name,
        ),
    )
    slurm_resolution_meta = {}
    if resolution is not None:
        slurm_resolution_meta = {
            "requested_queue_name": resolution.requested_queue,
            "resolved_queue_name": resolution.resolved_queue,
            "used_fallback": resolution.used_fallback,
            "slurm_partition": resolution.mapping.partition,
            "slurm_account": resolution.mapping.account,
            "slurm_qos": resolution.mapping.qos,
            "slurm_max_walltime_minutes": resolution.mapping.max_walltime_minutes,
        }
    meta_store = get_job_meta_store()
    meta_store.set_meta(
        job_id,
        {
            "request_id": request_id,
            "user_id": user_id,
            "tenant_id": tenant_id,
            "structure_id": request.structure_id,
            "queue_name": resolved_queue_name,
            "calc_mode": request.calc_mode,
            **slurm_resolution_meta,
        },
    )
    return SubmitJobOutcome(
        job_id=job_id,
        requested_queue_name=target.queue_name,
        resolved_queue_name=resolved_queue_name,
        slurm_resolution=resolution,
        idempotent=False,
    )


def get_job_status(job_id: str) -> ZPEJobStatus:
    store = get_result_store()
    return store.get_status(job_id)


def _ensure_job_finished(store: ResultStore, job_id: str) -> None:
    try:
        status = store.get_status(job_id)
    except KeyError as exc:
        raise KeyError("job not found") from exc
    if status.status == "failed":
        detail = status.detail or "unknown"
        raise JobConflictError(f"job failed: {detail}")
    if status.status != "finished":
        raise JobConflictError(f"job not finished (status={status.status})")


def get_job_result(job_id: str) -> Dict[str, Any]:
    store = get_result_store()
    _ensure_job_finished(store, job_id)
    try:
        result = store.get_result(job_id)
    except KeyError as exc:
        raise ResultArtifactMissingError("result missing after completion") from exc
    if "kpts" in result and "kpoints" not in result:
        result = {
            **{key: value for key, value in result.items() if key != "kpts"},
            "kpoints": result.get("kpts"),
        }
    return result


def get_job_file(job_id: str, kind: str) -> str:
    store = get_result_store()
    _ensure_job_finished(store, job_id)
    try:
        return store.get_file(job_id, kind)
    except KeyError as exc:
        raise ResultArtifactMissingError("file missing after completion") from exc


def submit_compute_result(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    event_id: str | None,
    result: Dict[str, Any],
    summary_text: str,
    freqs_csv: str,
    worker_meta: Dict[str, Any] | None = None,
) -> SubmitComputeResultOutcome:
    job_meta = get_job_meta_store().get_meta(job_id)
    outcome = submit_result(
        job_id=job_id,
        worker_id=worker_id,
        lease_id=lease_id,
        event_id=event_id,
        result=result,
        summary_text=summary_text,
        freqs_csv=freqs_csv,
    )
    duration_ms: int | None = None
    qe_version = worker_meta.get("qe_version") if worker_meta else None
    if worker_meta and "computation_time_seconds" in worker_meta:
        try:
            duration_ms = int(float(worker_meta["computation_time_seconds"]) * 1000)
        except (TypeError, ValueError):
            duration_ms = None
    return SubmitComputeResultOutcome(
        outcome=outcome,
        job_meta=job_meta,
        duration_ms=duration_ms,
        qe_version=qe_version if isinstance(qe_version, str) else None,
    )


def submit_compute_failure(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    event_id: str | None,
    error_code: str,
    error_message: str,
    traceback: str | None = None,
) -> SubmitComputeFailureOutcome:
    job_meta = get_job_meta_store().get_meta(job_id)
    outcome = submit_failure(
        job_id=job_id,
        worker_id=worker_id,
        lease_id=lease_id,
        event_id=event_id,
        error_code=error_code,
        error_message=error_message,
        traceback=traceback,
    )
    return SubmitComputeFailureOutcome(outcome=outcome, job_meta=job_meta)


def enqueue_zpe_job(payload: Dict[str, Any], *, queue_name: str | None = None) -> str:
    settings = get_zpe_settings()
    if settings.compute_mode not in {"remote-queue", "remote-http", "mock"}:
        raise ValueError(
            "compute_mode must be 'remote-queue', 'remote-http', or 'mock'."
        )
    store = get_result_store()
    if settings.compute_mode == "remote-queue":
        job = get_queue(queue_name).enqueue(
            "services.zpe.worker.run_zpe_job",
            payload,
            job_timeout=settings.job_timeout_seconds,
            result_ttl=settings.result_ttl_seconds,
        )
        job_id = job.id
        if not isinstance(job_id, str) or not job_id:
            raise RuntimeError("RQ enqueue returned an invalid job id")
        store.set_status(job_id, "queued")
        return job_id
    if settings.compute_mode == "remote-http":
        return enqueue_http_job(payload)
    return _run_mock_job(payload, store)


def _run_mock_job(payload: Dict[str, Any], store: ResultStore) -> str:
    request = ZPEJobRequest(**payload)
    job_id = f"mock-{uuid4().hex}"
    store.set_status(job_id, "started")

    fixed_indices = extract_fixed_indices(request.content)
    kpts = parse_kpoints_automatic(request.content)
    kpts_use = kpts[0] if kpts else _default_kpts()

    nfreq = max(1, len(request.mobile_indices) * 3)
    freqs_cm = [100.0 + 5.0 * i for i in range(nfreq)]
    zpe_ev, s_vib_jmol_k = calc_zpe_and_s_vib(
        freqs_cm,
        low_cut_cm=get_zpe_settings().low_cut_cm,
        temperature=get_zpe_settings().temperature,
    )

    now = _now_iso()
    result: Dict[str, Any] = {
        "calc_type": request.calc_type,
        "freqs_cm": normalize_frequencies(freqs_cm),
        "zpe_ev": zpe_ev,
        "s_vib_jmol_k": s_vib_jmol_k,
        "mobile_indices": request.mobile_indices,
        "fixed_indices": fixed_indices,
        "kpts": kpts_use,
        "delta": get_zpe_settings().delta,
        "low_cut_cm": get_zpe_settings().low_cut_cm,
        "temperature": get_zpe_settings().temperature,
        "use_environ": request.use_environ,
        "qe_input": "mock",
        "pseudo_dir": "mock",
        "calc_start_time": now,
        "calc_end_time": now,
        "elapsed_seconds": 0.0,
        "cache_checked": 0,
        "cache_deleted": 0,
        "ecutwfc": None,
        "ecutrho": None,
    }

    summary_text = format_summary(
        result,
        pseudo_dir="mock",
        qe_input="mock",
        new_calc=request.calc_mode == "new",
    )
    freqs_csv = format_freqs_csv(freqs_cm)
    store.set_result(job_id, result, summary_text=summary_text, freqs_csv=freqs_csv)
    store.set_status(job_id, "finished")
    return job_id
