from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from app.schemas.zpe import ComputeFailedRequest, ComputeResultRequest


def _openapi_schema() -> dict[str, Any]:
    repo_root = Path(__file__).resolve().parents[3]
    schema_path = repo_root / "packages" / "api-client" / "openapi" / "openapi.json"
    return json.loads(schema_path.read_text(encoding="utf-8"))


def test_openapi_contract_covers_submit_running_terminal_payloads() -> None:
    schema = _openapi_schema()
    components = schema["components"]["schemas"]  # type: ignore[index]
    paths = schema["paths"]  # type: ignore[index]

    submit_command = components["SubmitJobCommand"]  # type: ignore[index]
    assert set(submit_command["required"]) >= {  # type: ignore[index]
        "tenant_id",
        "workspace_id",
        "job_id",
        "idempotency_key",
        "management_node_id",
        "execution_profile",
        "resource_shape",
        "payload_ref",
        "requested_by",
    }
    event_payload = components["ExecutionEvent"]  # type: ignore[index]
    assert set(event_payload["required"]) >= {  # type: ignore[index]
        "event_id",
        "tenant_id",
        "workspace_id",
        "job_id",
        "submission_id",
        "execution_id",
        "state",
        "occurred_at",
        "trace_id",
    }

    runtime_status = components["RuntimeJobStatusResponse"]  # type: ignore[index]
    assert set(runtime_status["required"]) >= {  # type: ignore[index]
        "job_id",
        "submission_id",
        "execution_owner",
        "state",
        "updated_at",
        "trace_id",
    }

    ack = components["RuntimeEventAck"]  # type: ignore[index]
    assert "ok" in ack["properties"]  # type: ignore[index]
    assert "idempotent" in ack["properties"]  # type: ignore[index]

    assert "/api/runtime/jobs:submit" in paths
    assert "/api/runtime/jobs/{job_id}/events" in paths
    assert "/api/runtime/jobs/{job_id}" in paths
    assert "/api/runtime/jobs/{job_id}/projection" in paths


def test_management_node_payload_examples_validate_against_contract_models() -> None:
    result_payload: dict[str, Any] = {
        "tenant_id": "tenant-contract",
        "lease_id": "lease-123",
        "result": {"calc_type": "qe.zpe.v1", "zpe_ev": 0.12},
        "summary_text": "ok",
        "freqs_csv": "frequency_cm^-1,intensity\\n100,1.0",
        "meta": {"worker_hostname": "node-1"},
    }
    failed_payload: dict[str, Any] = {
        "tenant_id": "tenant-contract",
        "lease_id": "lease-123",
        "error_code": "ERR_RUNTIME",
        "error_message": "boom",
        "traceback": "trace",
    }

    result = ComputeResultRequest(**result_payload)
    failed = ComputeFailedRequest(**failed_payload)

    assert result.tenant_id == "tenant-contract"
    assert failed.error_code == "ERR_RUNTIME"


def test_result_request_normalizes_aiida_runtime_execution_event_shape() -> None:
    payload: dict[str, Any] = {
        "tenant_id": "tenant-contract",
        "lease_id": "lease-123",
        "result": {"calc_type": "qe.zpe.v1", "zpe_ev": 0.12},
        "summary_text": "ok",
        "freqs_csv": "frequency_cm^-1,intensity\\n100,1.0",
        "execution_event": {
            "eventId": "evt-1",
            "tenantId": "tenant-contract",
            "workspaceId": "workspace-1",
            "jobId": "job-1",
            "submissionId": "submission-1",
            "executionId": "exec-1",
            "status": "SUCCEEDED",
            "timestamp": "2026-01-01T00:00:00+00:00",
            "traceId": "trace-1",
            "statusDetail": "done",
            "schedulerRef": {
                "slurmJobId": "12345",
                "partition": "short",
                "qos": "normal",
            },
            "resultRef": {
                "outputUri": "zpe://jobs/job-1/result",
                "metadataUri": "zpe://jobs/job-1/files/summary",
            },
        },
    }

    normalized = ComputeResultRequest(**payload)
    assert normalized.execution_event is not None
    assert normalized.execution_event.event_id == "evt-1"
    assert normalized.execution_event.state == "completed"
    assert normalized.execution_event.occurred_at == "2026-01-01T00:00:00+00:00"
    assert normalized.execution_event.trace_id == "trace-1"
    assert normalized.execution_event.scheduler_ref is not None
    assert normalized.execution_event.scheduler_ref.slurm_job_id == "12345"
    assert normalized.execution_event.result_ref.output_uri == "zpe://jobs/job-1/result"


def test_failed_request_normalizes_aiida_runtime_execution_event_shape() -> None:
    payload: dict[str, Any] = {
        "tenant_id": "tenant-contract",
        "lease_id": "lease-123",
        "error_code": "COMPUTE_ERROR",
        "error_message": "boom",
        "execution_event": {
            "eventId": "evt-2",
            "tenantId": "tenant-contract",
            "workspaceId": "workspace-1",
            "jobId": "job-1",
            "submissionId": "submission-1",
            "executionId": "exec-1",
            "status": "ERROR",
            "occurredAt": "2026-01-01T00:00:01+00:00",
            "traceId": "trace-2",
            "errorInfo": {
                "errorCode": "COMPUTE_ERROR",
                "errorMessage": "boom",
                "retriable": True,
            },
        },
    }

    normalized = ComputeFailedRequest(**payload)
    assert normalized.execution_event is not None
    assert normalized.execution_event.event_id == "evt-2"
    assert normalized.execution_event.state == "failed"
    assert normalized.execution_event.occurred_at == "2026-01-01T00:00:01+00:00"
    assert normalized.execution_event.error.code == "COMPUTE_ERROR"
    assert normalized.execution_event.error.message == "boom"
    assert normalized.execution_event.error.retryable is True


def test_failed_request_accepts_top_level_error_aliases() -> None:
    payload: dict[str, Any] = {
        "tenant_id": "tenant-contract",
        "lease_id": "lease-123",
        "errorCode": "COMPUTE_ERROR",
        "errorMessage": "boom",
        "retryable": True,
        "execution_event": {
            "eventId": "evt-3",
            "tenantId": "tenant-contract",
            "workspaceId": "workspace-1",
            "jobId": "job-1",
            "submissionId": "submission-1",
            "executionId": "exec-1",
            "status": "ERROR",
            "occurredAt": "2026-01-01T00:00:02+00:00",
            "traceId": "trace-3",
        },
    }

    normalized = ComputeFailedRequest(**payload)
    assert normalized.error_code == "COMPUTE_ERROR"
    assert normalized.error_message == "boom"
    assert normalized.execution_event is not None
    assert normalized.execution_event.error.code == "COMPUTE_ERROR"
    assert normalized.execution_event.error.message == "boom"
    assert normalized.execution_event.error.retryable is True
