from __future__ import annotations

from fastapi.testclient import TestClient

from app.api import app
from app.schemas.runtime import RuntimeJobStatusResponse, SubmitJobAccepted
from app.schemas.zpe import ZPEJobStatus
from services.runtime import (
    RuntimeConfigurationError,
    RuntimeConflictError,
    RuntimeDownstreamError,
    RuntimeNotFoundError,
)


def _headers(*, user_id: str = "user-1", tenant_id: str = "tenant-1") -> dict[str, str]:
    return {
        "Authorization": f"Bearer {user_id}",
        "X-Tenant-Id": tenant_id,
    }


def _submit_payload() -> dict[str, object]:
    return {
        "tenant_id": "tenant-1",
        "workspace_id": "ws-1",
        "job_id": "job-1",
        "idempotency_key": "idem-1",
        "management_node_id": "mn-1",
        "execution_profile": {"queue_name": "standard"},
        "resource_shape": {"cpu": 2, "memory_mib": 1024, "walltime_seconds": 3600},
        "payload_ref": {"input_uri": "s3://bucket/input.in"},
        "requested_by": {"user_id": "user-1"},
    }


class _GatewaySuccess:
    def submit(self, _command, *, trace_id):  # type: ignore[no-untyped-def]
        return type(
            "SubmitResult",
            (),
            {
                "response": SubmitJobAccepted(
                    job_id="job-1",
                    submission_id="sub-1",
                    execution_owner="aiida-slurm",
                    accepted_at="2026-02-18T00:00:00+00:00",
                    trace_id=trace_id,
                ),
                "idempotent": False,
            },
        )()

    def apply_event(self, _event):  # type: ignore[no-untyped-def]
        return type("EventAck", (), {"idempotent": False})()

    def get_status(self, _job_id, *, tenant_id):  # type: ignore[no-untyped-def]
        return RuntimeJobStatusResponse(
            job_id="job-1",
            submission_id="sub-1",
            execution_owner="aiida-slurm",
            state="running",
            detail=f"tenant={tenant_id}",
            updated_at="2026-02-18T00:00:00+00:00",
            trace_id="trace-1",
        )

    def get_detail(self, _job_id, *, tenant_id):  # type: ignore[no-untyped-def]
        return {"tenant_id": tenant_id, "state": "running"}

    def get_projection_status(self, _job_id, *, tenant_id):  # type: ignore[no-untyped-def]
        return ZPEJobStatus(
            status="started",
            detail=f"tenant={tenant_id}",
            updated_at="2026-02-18T00:00:00+00:00",
        )


def test_runtime_submit_success(monkeypatch):
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewaySuccess())
    client = TestClient(app)

    response = client.post("/api/runtime/jobs:submit", json=_submit_payload(), headers=_headers())

    assert response.status_code == 202
    assert response.json()["job_id"] == "job-1"
    assert response.json()["submission_id"] == "sub-1"


def test_runtime_submit_scope_violation(monkeypatch):
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewaySuccess())
    client = TestClient(app)
    payload = _submit_payload()
    payload["tenant_id"] = "tenant-other"

    response = client.post("/api/runtime/jobs:submit", json=payload, headers=_headers())
    assert response.status_code == 403


def test_runtime_submit_conflict_and_config_errors(monkeypatch):
    class _GatewayConflict(_GatewaySuccess):
        def submit(self, _command, *, trace_id):  # type: ignore[no-untyped-def]
            _ = trace_id
            raise RuntimeConflictError("idempotency_conflict")

    class _GatewayConfig(_GatewaySuccess):
        def submit(self, _command, *, trace_id):  # type: ignore[no-untyped-def]
            _ = trace_id
            raise RuntimeConfigurationError("missing runtime setting: RUNTIME_COMMAND_SUBMIT_URL")

    client = TestClient(app)

    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewayConflict())
    conflict = client.post("/api/runtime/jobs:submit", json=_submit_payload(), headers=_headers())
    assert conflict.status_code == 409

    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewayConfig())
    config = client.post("/api/runtime/jobs:submit", json=_submit_payload(), headers=_headers())
    assert config.status_code == 503


def test_runtime_event_and_reads(monkeypatch):
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewaySuccess())
    client = TestClient(app)

    event_payload = {
        "event_id": "evt-1",
        "tenant_id": "tenant-1",
        "workspace_id": "ws-1",
        "job_id": "job-1",
        "submission_id": "sub-1",
        "execution_id": "exec-1",
        "state": "running",
        "occurred_at": "2026-02-17T00:00:00+00:00",
        "trace_id": "trace-1",
    }

    event = client.post("/api/runtime/jobs/job-1/events", json=event_payload, headers=_headers())
    assert event.status_code == 200
    assert event.json()["idempotent"] is False

    status = client.get("/api/runtime/jobs/job-1", headers=_headers())
    assert status.status_code == 200
    assert status.json()["state"] == "running"

    detail = client.get("/api/runtime/jobs/job-1/detail", headers=_headers())
    assert detail.status_code == 200
    assert detail.json()["tenant_id"] == "tenant-1"

    projection = client.get("/api/runtime/jobs/job-1/projection", headers=_headers())
    assert projection.status_code == 200
    assert projection.json()["status"] == "started"


def test_runtime_reads_map_errors(monkeypatch):
    class _GatewayMissing(_GatewaySuccess):
        def get_status(self, _job_id, *, tenant_id):  # type: ignore[no-untyped-def]
            _ = tenant_id
            raise RuntimeNotFoundError("missing")

    class _GatewayDownstream(_GatewaySuccess):
        def get_status(self, _job_id, *, tenant_id):  # type: ignore[no-untyped-def]
            _ = tenant_id
            raise RuntimeDownstreamError("runtime downstream request failed")

    client = TestClient(app)

    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewayMissing())
    missing = client.get("/api/runtime/jobs/job-1", headers=_headers())
    assert missing.status_code == 404

    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewayDownstream())
    downstream = client.get("/api/runtime/jobs/job-1", headers=_headers())
    assert downstream.status_code == 502


def test_runtime_event_maps_configuration_error(monkeypatch):
    class _GatewayConfig(_GatewaySuccess):
        def apply_event(self, _event):  # type: ignore[no-untyped-def]
            raise RuntimeConfigurationError(
                "missing runtime setting: RUNTIME_COMMAND_EVENT_URL_TEMPLATE"
            )

    client = TestClient(app)
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: _GatewayConfig())
    payload = {
        "event_id": "evt-1",
        "tenant_id": "tenant-1",
        "workspace_id": "ws-1",
        "job_id": "job-1",
        "submission_id": "sub-1",
        "execution_id": "exec-1",
        "state": "running",
        "occurred_at": "2026-02-17T00:00:00+00:00",
        "trace_id": "trace-1",
    }
    response = client.post("/api/runtime/jobs/job-1/events", json=payload, headers=_headers())
    assert response.status_code == 503
