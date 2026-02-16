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
