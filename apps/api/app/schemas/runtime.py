from __future__ import annotations

from typing import Any, Literal
from typing import cast

from pydantic import Field, model_validator

from .base import ApiModel
from .zpe import _normalize_aiida_state, _normalize_execution_event_payload


RuntimeState = Literal["accepted", "running", "completed", "failed"]


class SubmitExecutionProfile(ApiModel):
    queue_name: str
    qos: str | None = None
    account: str | None = None


class SubmitResourceShape(ApiModel):
    cpu: int = Field(ge=1)
    memory_mib: int = Field(ge=1)
    walltime_seconds: int = Field(ge=1)


class SubmitPayloadRef(ApiModel):
    input_uri: str
    artifact_bucket: str | None = None


class SubmitRequestedBy(ApiModel):
    user_id: str
    session_id: str | None = None


class SubmitJobCommand(ApiModel):
    tenant_id: str
    workspace_id: str
    job_id: str
    idempotency_key: str
    management_node_id: str
    execution_profile: SubmitExecutionProfile
    resource_shape: SubmitResourceShape
    payload_ref: SubmitPayloadRef
    requested_by: SubmitRequestedBy


class SubmitJobAccepted(ApiModel):
    job_id: str
    submission_id: str
    execution_owner: str
    accepted_at: str
    trace_id: str


class SchedulerRef(ApiModel):
    slurm_job_id: str | None = None
    partition: str | None = None
    qos: str | None = None


class ResultRef(ApiModel):
    output_uri: str
    metadata_uri: str | None = None


class EventError(ApiModel):
    code: str
    message: str
    retryable: bool


class ExecutionEvent(ApiModel):
    event_id: str
    tenant_id: str
    workspace_id: str
    job_id: str
    submission_id: str
    execution_id: str
    state: RuntimeState
    occurred_at: str
    trace_id: str
    status_detail: str | None = None
    scheduler_ref: SchedulerRef | None = None
    result_ref: ResultRef | None = None
    error: EventError | None = None

    @model_validator(mode="before")
    @classmethod
    def _normalize_event(cls, value: Any) -> Any:
        payload = _normalize_execution_event_payload(value)
        if isinstance(payload, dict) and "state" in payload:
            normalized = cast(dict[str, Any], payload)
            normalized["state"] = _normalize_aiida_state(normalized.get("state"))
            return normalized
        return cast(Any, payload)


class RuntimeJobStatusResponse(ApiModel):
    job_id: str
    submission_id: str
    execution_owner: str
    state: RuntimeState
    detail: str | None = None
    updated_at: str
    trace_id: str


class RuntimeEventAck(ApiModel):
    ok: bool = True
    idempotent: bool = False
