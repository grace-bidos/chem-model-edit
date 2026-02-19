from __future__ import annotations

from types import SimpleNamespace

from fastapi.testclient import TestClient
import pytest

from app.api import app
from app.schemas.runtime import ExecutionEvent, RuntimeJobStatusResponse, SubmitJobAccepted
from app.schemas.zpe import ZPEJobStatus
from services.runtime import (
    RuntimeGateway,
    RuntimeConfigurationError,
    RuntimeConflictError,
    RuntimeDownstreamError,
    RuntimeNotFoundError,
)
from services.runtime_nodes import RuntimeJoinToken, RuntimeTarget


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


def test_runtime_gateway_apply_event_posts_to_command_endpoint(monkeypatch):
    settings = SimpleNamespace(
        request_timeout_seconds=5,
        service_auth_bearer_token=None,
        command_submit_url="https://gateway.example/api/runtime/jobs:submit",
        command_event_url_template="https://gateway.example/api/runtime/jobs/{job_id}/events",
        read_status_url_template="https://gateway.example/api/runtime/jobs/{job_id}",
        read_detail_url_template="https://gateway.example/api/runtime/jobs/{job_id}/detail",
        read_projection_url_template="https://gateway.example/api/runtime/jobs/{job_id}/projection",
    )
    monkeypatch.setattr("services.runtime.get_runtime_settings", lambda: settings)
    gateway = RuntimeGateway()
    captured: dict[str, object] = {}

    def _request_json(**kwargs):  # type: ignore[no-untyped-def]
        captured.update(kwargs)
        return {"idempotent": True}

    monkeypatch.setattr(gateway, "_request_json", _request_json)
    event = ExecutionEvent(
        event_id="evt-1",
        tenant_id="tenant-1",
        workspace_id="ws-1",
        job_id="job/1",
        submission_id="sub-1",
        execution_id="exec-1",
        state="running",
        occurred_at="2026-02-17T00:00:00+00:00",
        trace_id="trace-1",
    )
    ack = gateway.apply_event(event)

    assert ack.idempotent is True
    assert captured["method"] == "POST"
    assert captured["url"] == "https://gateway.example/api/runtime/jobs/job%2F1/events"
    assert captured["tenant_id"] == "tenant-1"


def test_runtime_gateway_apply_event_requires_event_template(monkeypatch):
    settings = SimpleNamespace(
        request_timeout_seconds=5,
        service_auth_bearer_token=None,
        command_submit_url="https://gateway.example/api/runtime/jobs:submit",
        command_event_url_template=None,
        read_status_url_template="https://gateway.example/api/runtime/jobs/{job_id}",
        read_detail_url_template="https://gateway.example/api/runtime/jobs/{job_id}/detail",
        read_projection_url_template="https://gateway.example/api/runtime/jobs/{job_id}/projection",
    )
    monkeypatch.setattr("services.runtime.get_runtime_settings", lambda: settings)
    gateway = RuntimeGateway()
    event = ExecutionEvent(
        event_id="evt-1",
        tenant_id="tenant-1",
        workspace_id="ws-1",
        job_id="job-1",
        submission_id="sub-1",
        execution_id="exec-1",
        state="running",
        occurred_at="2026-02-17T00:00:00+00:00",
        trace_id="trace-1",
    )

    with pytest.raises(RuntimeConfigurationError, match="RUNTIME_COMMAND_EVENT_URL_TEMPLATE"):
        gateway.apply_event(event)


def test_runtime_issue_join_token(monkeypatch):
    class _RuntimeNodeStore:
        def create_join_token(  # type: ignore[no-untyped-def]
            self, *, tenant_id, owner_user_id, queue_name, ttl_seconds, node_name_hint
        ):
            _ = tenant_id
            _ = owner_user_id
            _ = ttl_seconds
            return RuntimeJoinToken(
                token="join-token-1",
                tenant_id="tenant-1",
                owner_user_id="user-1",
                queue_name=queue_name,
                created_at="2026-02-18T00:00:00+00:00",
                expires_at="2026-02-18T00:30:00+00:00",
                ttl_seconds=1800,
                node_name_hint=node_name_hint,
                label=queue_name,
            )

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: _RuntimeNodeStore())
    client = TestClient(app)
    response = client.post(
        "/api/runtime/nodes/join-token",
        json={"queue_name": "highmem", "ttl_seconds": 1800, "node_name_hint": "node-1"},
        headers=_headers(),
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["join_token"] == "join-token-1"
    assert payload["queue_name"] == "highmem"
    assert "--join-token join-token-1" in payload["install_command"]
    assert "--queue-name highmem" in payload["install_command"]


def test_runtime_install_script_endpoint_available():
    client = TestClient(app)
    response = client.get("/api/runtime/nodes/install.sh")
    assert response.status_code == 200
    assert "install-compute-node.sh" in response.text
    assert "--join-token" in response.text


def test_runtime_register_compute_node(monkeypatch):
    class _RuntimeNodeStore:
        def __init__(self):
            self.active_target: RuntimeTarget | None = None

        def consume_join_token(self, _token):  # type: ignore[no-untyped-def]
            return RuntimeJoinToken(
                token="join-token-1",
                tenant_id="tenant-1",
                owner_user_id="user-1",
                queue_name="standard",
                created_at="2026-02-18T00:00:00+00:00",
                expires_at="2026-02-18T00:30:00+00:00",
                ttl_seconds=1800,
                label="standard",
            )

        def add_target(self, *, tenant_id, user_id, queue_name, server_id, name=None, metadata=None):  # type: ignore[no-untyped-def]
            _ = tenant_id
            _ = user_id
            _ = server_id
            _ = metadata
            return RuntimeTarget(
                target_id="qt-1",
                tenant_id="tenant-1",
                user_id="user-1",
                queue_name=queue_name,
                server_id="compute-1",
                registered_at="2026-02-18T00:45:00+00:00",
                name=name,
            )

        def get_active_target(self, _tenant_id, _user_id):  # type: ignore[no-untyped-def]
            return self.active_target

        def set_active_target(self, _tenant_id, _user_id, _target_id):  # type: ignore[no-untyped-def]
            self.active_target = RuntimeTarget(
                target_id="qt-1",
                tenant_id="tenant-1",
                user_id="user-1",
                queue_name="standard",
                server_id="compute-1",
                registered_at="2026-02-18T00:45:00+00:00",
            )

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: _RuntimeNodeStore())
    client = TestClient(app)
    response = client.post(
        "/api/runtime/nodes/register",
        json={"token": "join-token-1", "name": "worker-1"},
    )
    assert response.status_code == 200
    payload = response.json()
    assert payload["server_id"].startswith("compute-")
    assert payload["target_id"] == "qt-1"
    assert payload["queue_name"] == "standard"


def test_runtime_register_requires_owner(monkeypatch):
    class _RuntimeNodeStore:
        def consume_join_token(self, _token):  # type: ignore[no-untyped-def]
            raise KeyError("join token not found")

    monkeypatch.setattr("app.routers.runtime.get_runtime_node_store", lambda: _RuntimeNodeStore())
    client = TestClient(app)
    response = client.post(
        "/api/runtime/nodes/register",
        json={"token": "join-token-1"},
    )
    assert response.status_code == 404
