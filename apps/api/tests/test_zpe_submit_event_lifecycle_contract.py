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

    lease_response = components["ComputeLeaseResponse"]  # type: ignore[index]
    assert set(lease_response["required"]) >= {  # type: ignore[index]
        "job_id",
        "payload",
        "lease_id",
        "lease_ttl_seconds",
    }
    assert "meta" in lease_response["properties"]  # type: ignore[index]

    result_request = components["ComputeResultRequest"]  # type: ignore[index]
    assert set(result_request["required"]) >= {  # type: ignore[index]
        "tenant_id",
        "lease_id",
        "result",
        "summary_text",
        "freqs_csv",
    }

    failed_request = components["ComputeFailedRequest"]  # type: ignore[index]
    assert set(failed_request["required"]) >= {  # type: ignore[index]
        "tenant_id",
        "lease_id",
        "error_code",
        "error_message",
    }

    status_schema = components["ZPEJobStatus"]  # type: ignore[index]
    state_enum = status_schema["properties"]["status"]["enum"]  # type: ignore[index]
    assert {"queued", "started", "finished", "failed"} <= set(state_enum)

    assert "/api/zpe/jobs" in paths
    assert "/api/zpe/compute/jobs/lease" in paths
    assert "/api/zpe/compute/jobs/{job_id}/result" in paths
    assert "/api/zpe/compute/jobs/{job_id}/failed" in paths


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
