from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Tuple, cast

from pydantic import Field, model_validator

from .base import ApiModel
from .common import Pagination, Structure


JobState = Literal["queued", "started", "finished", "failed"]


_AIIA_STATE_MAP: dict[str, str] = {
    "complete": "completed",
    "completed": "completed",
    "done": "completed",
    "error": "failed",
    "failed": "failed",
    "failure": "failed",
    "finished": "completed",
    "success": "completed",
    "succeeded": "completed",
}


def _to_string_key_dict(value: Any) -> dict[str, Any] | None:
    if not isinstance(value, dict):
        return None
    source = cast(dict[object, Any], value)
    normalized: dict[str, Any] = {}
    for raw_key, raw_item in source.items():
        if isinstance(raw_key, str):
            normalized[raw_key] = raw_item
    return normalized


def _first_present(payload: dict[str, Any], *keys: str) -> Any:
    for key in keys:
        value = payload.get(key)
        if value is not None:
            return value
    return None


def _normalize_aiida_state(value: Any) -> Any:
    if not isinstance(value, str):
        return value
    normalized = value.strip().lower()
    return _AIIA_STATE_MAP.get(normalized, normalized)


def _normalize_scheduler_ref_payload(value: Any) -> Any:
    payload = _to_string_key_dict(value)
    if payload is None:
        return value
    slurm_job_id = _first_present(payload, "slurm_job_id", "slurmJobId")
    partition = _first_present(payload, "partition")
    qos = _first_present(payload, "qos")
    normalized: dict[str, Any] = {}
    if slurm_job_id is not None:
        normalized["slurm_job_id"] = slurm_job_id
    if partition is not None:
        normalized["partition"] = partition
    if qos is not None:
        normalized["qos"] = qos
    return normalized or payload


def _normalize_result_ref_payload(value: Any) -> Any:
    payload = _to_string_key_dict(value)
    if payload is None:
        return value
    output_uri = _first_present(payload, "output_uri", "outputUri")
    metadata_uri = _first_present(payload, "metadata_uri", "metadataUri")
    normalized: dict[str, Any] = {}
    if output_uri is not None:
        normalized["output_uri"] = output_uri
    if metadata_uri is not None:
        normalized["metadata_uri"] = metadata_uri
    return normalized or payload


def _normalize_error_payload(value: Any) -> Any:
    payload = _to_string_key_dict(value)
    if payload is None:
        return value
    code = _first_present(payload, "code", "error_code", "errorCode")
    message = _first_present(payload, "message", "error_message", "errorMessage")
    retryable = _first_present(payload, "retryable", "retriable")
    normalized: dict[str, Any] = {}
    if code is not None:
        normalized["code"] = code
    if message is not None:
        normalized["message"] = message
    if retryable is not None:
        normalized["retryable"] = retryable
    return normalized or payload


def _normalize_execution_event_payload(value: Any) -> Any:
    payload = _to_string_key_dict(value)
    if payload is None:
        return value
    normalized: dict[str, Any] = {}
    fields: tuple[tuple[str, tuple[str, ...]], ...] = (
        ("event_id", ("event_id", "eventId")),
        ("tenant_id", ("tenant_id", "tenantId")),
        ("workspace_id", ("workspace_id", "workspaceId")),
        ("job_id", ("job_id", "jobId")),
        ("submission_id", ("submission_id", "submissionId")),
        ("execution_id", ("execution_id", "executionId")),
        ("occurred_at", ("occurred_at", "occurredAt", "timestamp", "event_time")),
        ("trace_id", ("trace_id", "traceId")),
        ("status_detail", ("status_detail", "statusDetail", "detail")),
    )
    for target, keys in fields:
        normalized_value = _first_present(payload, *keys)
        if normalized_value is not None:
            normalized[target] = normalized_value

    state = _first_present(payload, "state", "status", "lifecycle_state")
    if state is not None:
        normalized["state"] = _normalize_aiida_state(state)

    scheduler_ref = _first_present(payload, "scheduler_ref", "schedulerRef", "scheduler")
    if scheduler_ref is not None:
        normalized["scheduler_ref"] = _normalize_scheduler_ref_payload(scheduler_ref)
    result_ref = _first_present(payload, "result_ref", "resultRef")
    if result_ref is not None:
        normalized["result_ref"] = _normalize_result_ref_payload(result_ref)
    elif "output_uri" in payload or "outputUri" in payload:
        normalized["result_ref"] = _normalize_result_ref_payload(payload)

    error = _first_present(payload, "error", "error_info", "errorInfo")
    if error is not None:
        normalized["error"] = _normalize_error_payload(error)

    return normalized or payload


class ZPEParseRequest(ApiModel):
    content: str
    structure_id: Optional[str] = None


class ZPEParseResponse(ApiModel):
    structure: Structure
    fixed_indices: List[int]
    atomic_species: Dict[str, str] = Field(default_factory=dict)
    kpoints: Optional[Tuple[int, int, int]] = None


class ZPEJobRequest(ApiModel):
    calc_type: Literal["qe.zpe.v1", "qe.relax.v1"] = "qe.zpe.v1"
    content: str
    mobile_indices: List[int]
    use_environ: bool = False
    input_dir: Optional[str] = None
    calc_mode: Literal["new", "continue"] = "continue"
    structure_id: Optional[str] = None


class ZPEJobResponse(ApiModel):
    id: str


class ZPEJobStatus(ApiModel):
    status: JobState
    detail: Optional[str] = None
    updated_at: Optional[str] = None


class ZPEResult(ApiModel):
    calc_type: Literal["qe.zpe.v1", "qe.relax.v1"] = "qe.zpe.v1"
    freqs_cm: List[float]
    zpe_ev: float
    s_vib_jmol_k: float
    mobile_indices: List[int]
    fixed_indices: List[int]
    kpoints: Tuple[int, int, int]
    delta: float
    low_cut_cm: float
    temperature: float
    use_environ: bool
    qe_input: str
    pseudo_dir: str
    calc_start_time: str
    calc_end_time: str
    elapsed_seconds: float
    cache_checked: int
    cache_deleted: int
    ecutwfc: Optional[float] = None
    ecutrho: Optional[float] = None


class ZPEJobResultResponse(ApiModel):
    result: ZPEResult


class EnrollTokenRequest(ApiModel):
    ttl_seconds: Optional[int] = None
    label: Optional[str] = None


class EnrollTokenResponse(ApiModel):
    token: str
    expires_at: str
    ttl_seconds: int
    label: Optional[str] = None


class ComputeRegisterRequest(ApiModel):
    token: str
    name: Optional[str] = None
    queue_name: Optional[str] = None
    meta: Dict[str, Any] = Field(default_factory=dict)
    activate_target: bool = Field(default=False)


class ComputeRegisterResponse(ApiModel):
    id: str
    registered_at: str
    name: Optional[str] = None
    worker_token: str
    token_expires_at: str
    token_ttl_seconds: int


class ComputeRevokeResponse(ApiModel):
    revoked_count: int


class QueueTarget(ApiModel):
    id: str
    queue_name: str
    server_id: str
    registered_at: str
    name: Optional[str] = None


class QueueTargetListResponse(ApiModel):
    targets: List[QueueTarget]
    active_target_id: Optional[str] = None
    pagination: Pagination


class QueueTargetSelectResponse(ApiModel):
    active_target_id: str


class ComputeLeaseResponse(ApiModel):
    job_id: str
    payload: Dict[str, Any]
    lease_id: str
    lease_ttl_seconds: int
    meta: Dict[str, Any] = Field(default_factory=dict)


class ExecutionEventSchedulerRef(ApiModel):
    slurm_job_id: Optional[str] = None
    partition: Optional[str] = None
    qos: Optional[str] = None

    @model_validator(mode="before")
    @classmethod
    def _normalize_scheduler_ref(cls, value: Any) -> Any:
        return _normalize_scheduler_ref_payload(value)


class ExecutionEventResultRef(ApiModel):
    output_uri: str
    metadata_uri: Optional[str] = None

    @model_validator(mode="before")
    @classmethod
    def _normalize_result_ref(cls, value: Any) -> Any:
        return _normalize_result_ref_payload(value)


class ExecutionEventError(ApiModel):
    code: str
    message: str
    retryable: bool

    @model_validator(mode="before")
    @classmethod
    def _normalize_error(cls, value: Any) -> Any:
        return _normalize_error_payload(value)


class ExecutionEventBase(ApiModel):
    event_id: str
    tenant_id: str
    workspace_id: str
    job_id: str
    submission_id: str
    execution_id: str
    occurred_at: str
    trace_id: str
    status_detail: Optional[str] = None
    scheduler_ref: Optional[ExecutionEventSchedulerRef] = None

    @model_validator(mode="before")
    @classmethod
    def _normalize_event(cls, value: Any) -> Any:
        return _normalize_execution_event_payload(value)


class ExecutionCompletedEvent(ExecutionEventBase):
    state: Literal["completed"]
    result_ref: ExecutionEventResultRef


class ExecutionFailedEvent(ExecutionEventBase):
    state: Literal["failed"]
    error: ExecutionEventError


class ComputeResultRequest(ApiModel):
    tenant_id: str
    lease_id: str
    result: Dict[str, Any]
    summary_text: str
    freqs_csv: str
    execution_event: Optional[ExecutionCompletedEvent] = None
    meta: Dict[str, Any] = Field(default_factory=dict)

    @model_validator(mode="before")
    @classmethod
    def _normalize_result_execution_event(cls, value: Any) -> Any:
        payload = _to_string_key_dict(value)
        if payload is None:
            return value
        execution_event = payload.get("execution_event")
        if execution_event is not None:
            payload["execution_event"] = _normalize_execution_event_payload(
                execution_event
            )
        return payload


class ComputeResultResponse(ApiModel):
    ok: bool = True
    idempotent: bool = False


class ComputeFailedRequest(ApiModel):
    tenant_id: str
    lease_id: str
    error_code: str
    error_message: str
    traceback: Optional[str] = None
    execution_event: Optional[ExecutionFailedEvent] = None

    @model_validator(mode="before")
    @classmethod
    def _normalize_failed_execution_event(cls, value: Any) -> Any:
        payload = _to_string_key_dict(value)
        if payload is None:
            return value
        # Accept legacy top-level aliases while keeping ApiModel(extra="forbid") strict.
        if "error_code" not in payload and "errorCode" in payload:
            payload["error_code"] = payload["errorCode"]
        if "error_message" not in payload and "errorMessage" in payload:
            payload["error_message"] = payload["errorMessage"]
        payload.pop("errorCode", None)
        payload.pop("errorMessage", None)
        execution_event = payload.get("execution_event")
        retryable_alias = _first_present(payload, "retryable", "retriable")
        payload.pop("retryable", None)
        payload.pop("retriable", None)
        if execution_event is None:
            return payload
        normalized_event = _normalize_execution_event_payload(execution_event)
        normalized_event_payload = _to_string_key_dict(normalized_event)
        if (
            normalized_event_payload is not None
            and normalized_event_payload.get("error") is None
        ):
            error_code = _first_present(payload, "error_code", "errorCode")
            error_message = _first_present(payload, "error_message", "errorMessage")
            if error_code is not None and error_message is not None:
                error_payload: dict[str, Any] = {
                    "code": error_code,
                    "message": error_message,
                }
                if retryable_alias is not None:
                    error_payload["retryable"] = retryable_alias
                normalized_event_payload["error"] = error_payload
        payload["execution_event"] = normalized_event_payload or normalized_event
        return payload


class ComputeFailedResponse(ApiModel):
    ok: bool = True
    requeued: bool
    retry_count: int

