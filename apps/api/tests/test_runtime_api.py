from __future__ import annotations

from pathlib import Path

from fastapi.testclient import TestClient

from app.api import app
from services.runtime import RuntimeStore


def _headers(*, user_id: str = "user-1", tenant_id: str = "tenant-1") -> dict[str, str]:
    return {
        "X-Dev-User-Id": user_id,
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


def test_runtime_submit_is_idempotent(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: store)

    client = TestClient(app)
    payload = _submit_payload()

    first = client.post("/api/runtime/jobs:submit", json=payload, headers=_headers())
    second = client.post("/api/runtime/jobs:submit", json=payload, headers=_headers())

    assert first.status_code == 202
    assert second.status_code == 202
    assert first.json()["submission_id"] == second.json()["submission_id"]


def test_runtime_submit_scope_violation(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: store)

    client = TestClient(app)
    payload = _submit_payload()
    payload["tenant_id"] = "tenant-other"

    response = client.post("/api/runtime/jobs:submit", json=payload, headers=_headers())
    assert response.status_code == 403


def test_runtime_event_updates_status(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: store)

    client = TestClient(app)
    submit = client.post("/api/runtime/jobs:submit", json=_submit_payload(), headers=_headers())
    assert submit.status_code == 202
    submission_id = submit.json()["submission_id"]

    event_payload = {
        "event_id": "evt-1",
        "tenant_id": "tenant-1",
        "workspace_id": "ws-1",
        "job_id": "job-1",
        "submission_id": submission_id,
        "execution_id": "exec-1",
        "state": "running",
        "occurred_at": "2026-02-17T00:00:00+00:00",
        "trace_id": "trace-1",
    }

    event_response = client.post(
        "/api/runtime/jobs/job-1/events",
        json=event_payload,
        headers=_headers(),
    )
    assert event_response.status_code == 200

    status_response = client.get("/api/runtime/jobs/job-1", headers=_headers())
    assert status_response.status_code == 200
    assert status_response.json()["state"] == "running"

    projection_response = client.get("/api/runtime/jobs/job-1/projection", headers=_headers())
    assert projection_response.status_code == 200
    assert projection_response.json()["status"] == "started"


def test_runtime_status_isolated_by_tenant(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr("app.routers.runtime.get_runtime_store", lambda: store)

    client = TestClient(app)
    submit = client.post("/api/runtime/jobs:submit", json=_submit_payload(), headers=_headers())
    assert submit.status_code == 202

    response = client.get(
        "/api/runtime/jobs/job-1",
        headers=_headers(tenant_id="tenant-other"),
    )
    assert response.status_code == 404
