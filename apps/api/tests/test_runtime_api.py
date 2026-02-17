from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from fastapi.testclient import TestClient
import pytest

from app.api import app
from app.schemas.runtime import ExecutionEvent, SubmitJobCommand
from services.runtime import RuntimeStore
from services.runtime import RuntimeConflictError, RuntimeNotFoundError
from services.runtime import _event_time_iso, _to_job_state
from services.runtime_settings import resolve_runtime_db_path


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


def _submit_command(
    *,
    tenant_id: str = "tenant-1",
    workspace_id: str = "ws-1",
    job_id: str = "job-1",
    idempotency_key: str = "idem-1",
) -> SubmitJobCommand:
    return SubmitJobCommand.model_validate(
        {
            "tenant_id": tenant_id,
            "workspace_id": workspace_id,
            "job_id": job_id,
            "idempotency_key": idempotency_key,
            "management_node_id": "mn-1",
            "execution_profile": {"queue_name": "standard"},
            "resource_shape": {"cpu": 2, "memory_mib": 1024, "walltime_seconds": 3600},
            "payload_ref": {"input_uri": "s3://bucket/input.in"},
            "requested_by": {"user_id": "user-1"},
        }
    )


def _event(
    *,
    event_id: str = "evt-1",
    state: str = "running",
    tenant_id: str = "tenant-1",
    workspace_id: str = "ws-1",
    job_id: str = "job-1",
    submission_id: str = "sub-1",
    occurred_at: str = "2026-02-17T00:00:00+00:00",
) -> ExecutionEvent:
    return ExecutionEvent.model_validate(
        {
            "event_id": event_id,
            "tenant_id": tenant_id,
            "workspace_id": workspace_id,
            "job_id": job_id,
            "submission_id": submission_id,
            "execution_id": "exec-1",
            "state": state,
            "occurred_at": occurred_at,
            "trace_id": "trace-1",
            "status_detail": "ok",
        }
    )


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


def test_event_time_iso_normalizes_invalid_and_naive_values():
    invalid = _event_time_iso("not-a-date")
    assert invalid.endswith("+00:00")

    naive = _event_time_iso("2026-02-17T12:30:45")
    assert naive == "2026-02-17T12:30:45+00:00"


def test_to_job_state_covers_all_runtime_states():
    assert _to_job_state("accepted") == "queued"
    assert _to_job_state("running") == "started"
    assert _to_job_state("completed") == "finished"
    assert _to_job_state("failed") == "failed"


def test_runtime_store_submit_conflicts(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr(store, "_dispatch_projection", lambda _event: None)

    first = store.submit(_submit_command(), trace_id="trace-1")
    assert first.idempotent is False

    with pytest.raises(RuntimeConflictError, match="idempotency_conflict"):
        store.submit(
            _submit_command(idempotency_key="idem-1", job_id="job-other"),
            trace_id="trace-2",
        )

    with pytest.raises(RuntimeConflictError, match="job_id_conflict"):
        store.submit(
            _submit_command(idempotency_key="idem-2", job_id="job-1"),
            trace_id="trace-3",
        )


def test_runtime_store_apply_event_conflicts(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr(store, "_dispatch_projection", lambda _event: None)

    submitted = store.submit(_submit_command(), trace_id="trace-1")
    submission_id = submitted.response.submission_id

    with pytest.raises(RuntimeNotFoundError):
        store.apply_event(_event(job_id="job-missing", submission_id=submission_id))

    with pytest.raises(RuntimeConflictError, match="scope_conflict"):
        store.apply_event(_event(submission_id="sub-other"))

    first = store.apply_event(_event(event_id="evt-dup", submission_id=submission_id))
    assert first.idempotent is False
    replay = store.apply_event(_event(event_id="evt-dup", submission_id=submission_id))
    assert replay.idempotent is True

    with pytest.raises(RuntimeConflictError, match="event_id_conflict"):
        store.apply_event(
            _event(event_id="evt-dup", state="completed", submission_id=submission_id)
        )

    with pytest.raises(RuntimeConflictError, match="invalid_transition"):
        store.apply_event(
            _event(event_id="evt-backward", state="accepted", submission_id=submission_id)
        )

    complete = store.apply_event(
        _event(event_id="evt-complete", state="completed", submission_id=submission_id)
    )
    assert complete.idempotent is False

    with pytest.raises(RuntimeConflictError, match="invalid_transition"):
        store.apply_event(
            _event(event_id="evt-after-terminal", state="running", submission_id=submission_id)
        )


def test_runtime_store_get_detail_honors_tenant_scope(monkeypatch, tmp_path: Path):
    store = RuntimeStore(tmp_path / "runtime.sqlite3")
    monkeypatch.setattr(store, "_dispatch_projection", lambda _event: None)
    submitted = store.submit(_submit_command(), trace_id="trace-1")
    store.apply_event(_event(submission_id=submitted.response.submission_id))

    detail = store.get_detail("job-1", tenant_id="tenant-1")
    assert detail["last_event"] is not None

    with pytest.raises(RuntimeNotFoundError):
        store.get_detail("job-1", tenant_id="tenant-other")


def test_resolve_runtime_db_path_handles_absolute_and_relative(monkeypatch):
    monkeypatch.setattr(
        "services.runtime_settings.get_runtime_settings",
        lambda: SimpleNamespace(sqlite_path="/tmp/runtime.sqlite3"),
    )
    assert resolve_runtime_db_path() == Path("/tmp/runtime.sqlite3")

    monkeypatch.setattr(
        "services.runtime_settings.get_runtime_settings",
        lambda: SimpleNamespace(sqlite_path=".runtime/runtime.sqlite3"),
    )
    assert resolve_runtime_db_path() == Path.cwd() / ".runtime/runtime.sqlite3"
