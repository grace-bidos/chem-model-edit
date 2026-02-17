from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import json
import threading
import time
import fakeredis
from fastapi.testclient import TestClient
import pytest

import main
import app.routers.zpe as zpe_router
from app.schemas.zpe import ZPEJobRequest
from services.zpe import backends as zpe_backends
from services.zpe import compute_results as zpe_compute_results
from services.zpe import enroll as zpe_enroll
from services.zpe import job_meta as zpe_job_meta
from services.zpe import job_owner as zpe_job_owner
from services.zpe import ops_flags as zpe_ops_flags
from services.zpe import queue as zpe_queue
from services.zpe import queue_targets as zpe_queue_targets
from services.zpe import result_store as zpe_store
from services.zpe import submit_idempotency as zpe_submit_idempotency
from services.zpe import worker_auth as zpe_worker_auth
from services.zpe.result_store import RedisResultStore
from services.zpe.settings import ZPESettings

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=1, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
""".strip()


def _write_slurm_policy(path, *, fallback_mode: str = "route-default") -> None:
    policy = {
        "cluster_name": "chem-cluster",
        "slurm": {
            "partitions": ["short", "long"],
            "accounts": ["chem-default", "chem-premium"],
            "qos": ["normal", "priority"],
        },
        "queue_mappings": [
            {
                "queue": "standard",
                "partition": "short",
                "account": "chem-default",
                "qos": "normal",
                "max_walltime_minutes": 120,
            }
        ],
        "fallback_policy": {"mode": fallback_mode, "default_queue": "standard"},
    }
    if fallback_mode == "deny":
        policy["fallback_policy"] = {"mode": "deny"}
    path.write_text(json.dumps(policy), encoding="utf-8")


def _patch_redis(monkeypatch):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_store, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_compute_results, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_enroll, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_worker_auth, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_queue_targets, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_owner, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_job_meta, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_ops_flags, "get_redis_connection", lambda: fake)
    monkeypatch.setattr(zpe_submit_idempotency, "get_redis_connection", lambda: fake)
    return fake


def _setup_user_and_target(
    client: TestClient,
    monkeypatch,
    fake,
    *,
    user_id: str = "dev-user-1",
    queue_name: str = "zpe",
) -> tuple[dict[str, str], str]:
    tenant_id = f"tenant-{user_id}"
    headers = {
        "X-Dev-User-Id": user_id,
        "X-Dev-User-Email": f"{user_id}@example.com",
        "X-Tenant-Id": tenant_id,
    }

    enroll = client.post(
        "/api/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60},
        headers=headers,
    )
    enroll_token = enroll.json()["token"]

    client.post(
        "/api/zpe/compute/servers",
        json={
            "token": enroll_token,
            "name": "server-1",
            "queue_name": queue_name,
        },
    )
    return headers, user_id


def test_zpe_mock_api_flow(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    status = client.get(f"/api/zpe/jobs/{job_id}", headers=headers)
    assert status.status_code == 200
    assert status.json()["status"] == "finished"

    result = client.get(f"/api/zpe/jobs/{job_id}/result", headers=headers)
    assert result.status_code == 200
    payload = result.json()["result"]
    assert payload["zpe_ev"] >= 0.0
    assert payload["mobile_indices"] == [0]

    summary = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "summary"},
        headers=headers,
    )
    assert summary.status_code == 200
    assert "ZPE summary" in summary.text

    freqs = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "freqs"},
        headers=headers,
    )
    assert freqs.status_code == 200
    assert freqs.text.startswith("frequency_cm^-1")


def test_zpe_job_status_missing(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.get("/api/zpe/jobs/missing-job", headers=headers)
    assert response.status_code == 404


def test_zpe_job_result_not_finished(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    job_id = "job-not-finished"
    store.set_status(job_id, "queued")

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": headers["X-Tenant-Id"],
            "user_id": user_id,
            "request_id": "req-not-finished",
        },
    )
    response = client.get(f"/api/zpe/jobs/{job_id}/result", headers=headers)
    assert response.status_code == 409


def test_zpe_submission_disabled(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    zpe_ops_flags.set_ops_flags(submission_enabled=False, redis=fake)

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 503


def test_zpe_submission_route_next_gen_keeps_submit_hotpath_available(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    zpe_ops_flags.set_ops_flags(submission_route="next-gen", redis=fake)

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200


def test_zpe_submission_denies_unknown_queue_when_slurm_policy_is_deny(
    monkeypatch,
    tmp_path,
):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    policy_path = tmp_path / "onboarding.policy.json"
    _write_slurm_policy(policy_path, fallback_mode="deny")
    settings = ZPESettings(
        compute_mode="mock",
        result_store="redis",
        slurm_policy_path=str(policy_path),
    )

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(
        client,
        monkeypatch,
        fake,
        queue_name="legacy",
    )

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 400
    assert "not allowed by slurm fallback policy" in response.json()["error"]["message"]


def test_zpe_submission_surfaces_slurm_policy_read_error_as_503(
    monkeypatch,
    tmp_path,
):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    policy_path = tmp_path / "policy-dir"
    policy_path.mkdir()
    settings = ZPESettings(
        compute_mode="mock",
        result_store="redis",
        slurm_policy_path=str(policy_path),
    )

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 503
    assert "slurm policy configuration error" in response.json()["error"]["message"]


def test_zpe_submission_routes_to_default_queue_with_slurm_route_default_policy(
    monkeypatch,
    tmp_path,
):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    policy_path = tmp_path / "onboarding.policy.json"
    _write_slurm_policy(policy_path, fallback_mode="route-default")
    settings = ZPESettings(
        compute_mode="mock",
        result_store="redis",
        slurm_policy_path=str(policy_path),
    )
    captured: dict[str, str] = {}

    def _capture_enqueue(payload, *, queue_name=None):
        captured["queue_name"] = queue_name
        return "job-slurm-route-default"

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_backends, "enqueue_zpe_job", _capture_enqueue)

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(
        client,
        monkeypatch,
        fake,
        queue_name="legacy",
    )
    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200
    assert response.json()["id"] == "job-slurm-route-default"
    assert captured["queue_name"] == "standard"

    meta = zpe_job_meta.JobMetaStore(redis=fake).get_meta("job-slurm-route-default")
    assert meta["user_id"] == user_id
    assert meta["requested_queue_name"] == "legacy"
    assert meta["resolved_queue_name"] == "standard"
    assert meta["used_fallback"] is True
    assert meta["slurm_partition"] == "short"
    assert meta["slurm_account"] == "chem-default"
    assert meta["slurm_qos"] == "normal"
    assert meta["slurm_adapter_configured"] == "stub-policy"
    assert meta["slurm_adapter_effective"] == "stub-policy"
    assert meta["slurm_adapter_rollback_guard"] == "allow"


def test_zpe_submission_with_real_policy_surfaces_precondition_failure_as_503(
    monkeypatch,
    tmp_path,
):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    policy_path = tmp_path / "onboarding.policy.json"
    _write_slurm_policy(policy_path, fallback_mode="route-default")
    settings = ZPESettings(
        compute_mode="mock",
        result_store="redis",
        slurm_policy_path=str(policy_path),
        slurm_adapter="real-policy",
    )

    def _missing_scontrol(*args, **kwargs):
        _ = (args, kwargs)
        raise FileNotFoundError("scontrol")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr("services.zpe.slurm_policy.subprocess.run", _missing_scontrol)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 503
    assert "slurm real-adapter precondition failed" in response.json()["error"]["message"]


def test_zpe_submission_with_guard_forces_stub_policy_and_records_override(
    monkeypatch,
    tmp_path,
):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    policy_path = tmp_path / "onboarding.policy.json"
    _write_slurm_policy(policy_path, fallback_mode="route-default")
    settings = ZPESettings(
        compute_mode="mock",
        result_store="redis",
        slurm_policy_path=str(policy_path),
        slurm_adapter="real-policy",
        slurm_adapter_rollback_guard="force-stub-policy",
    )
    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(
        zpe_backends,
        "enqueue_zpe_job",
        lambda payload, *, queue_name=None: "job-guarded-real-policy",
    )
    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(
        client,
        monkeypatch,
        fake,
        queue_name="legacy",
    )

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200

    meta = zpe_job_meta.JobMetaStore(redis=fake).get_meta("job-guarded-real-policy")
    assert meta["slurm_adapter"] == "stub-policy"
    assert meta["slurm_adapter_configured"] == "real-policy"
    assert meta["slurm_adapter_effective"] == "stub-policy"
    assert meta["slurm_adapter_rollback_guard"] == "force-stub-policy"


def test_zpe_result_read_projection_keeps_read_hotpath_available(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    job_id = "job-projection-read"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": headers["X-Tenant-Id"],
            "user_id": user_id,
            "request_id": "req-projection-read",
        },
    )
    store.set_status(job_id, "finished")
    zpe_ops_flags.set_ops_flags(result_read_source="projection", redis=fake)

    status = client.get(f"/api/zpe/jobs/{job_id}", headers=headers)
    assert status.status_code == 200
    assert status.json()["status"] == "finished"


def test_admin_ops_flags_include_cutover_fields(monkeypatch, admin_auth_headers):
    _patch_redis(monkeypatch)
    client = TestClient(main.app)

    current = client.get("/api/zpe/admin/ops", headers=admin_auth_headers)
    assert current.status_code == 200
    payload = current.json()
    assert payload["submission_route"] == "redis-worker"
    assert payload["result_read_source"] == "redis"
    assert payload["legacy_worker_endpoints_enabled"] is True

    updated = client.patch(
        "/api/zpe/admin/ops",
        json={
            "submission_route": "next-gen",
            "result_read_source": "projection",
            "legacy_worker_endpoints_enabled": False,
        },
        headers=admin_auth_headers,
    )
    assert updated.status_code == 200
    out = updated.json()
    assert out["submission_route"] == "next-gen"
    assert out["result_read_source"] == "projection"
    assert out["legacy_worker_endpoints_enabled"] is False


def test_admin_ops_flags_partial_update_keeps_other_values(
    monkeypatch, admin_auth_headers
):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)

    first = client.patch(
        "/api/zpe/admin/ops",
        json={
            "submission_route": "next-gen",
            "result_read_source": "projection",
            "legacy_worker_endpoints_enabled": False,
        },
        headers=admin_auth_headers,
    )
    assert first.status_code == 200

    partial = client.patch(
        "/api/zpe/admin/ops",
        json={"dequeue_enabled": False},
        headers=admin_auth_headers,
    )
    assert partial.status_code == 200
    payload = partial.json()
    assert payload["dequeue_enabled"] is False
    assert payload["submission_route"] == "next-gen"
    assert payload["result_read_source"] == "projection"
    assert payload["legacy_worker_endpoints_enabled"] is False
    assert fake.ttl("zpe:ops:dequeue_enabled") > 0
    assert fake.ttl("zpe:ops:submission_route") > 0
    assert fake.ttl("zpe:ops:result_read_source") > 0
    assert fake.ttl("zpe:ops:legacy_worker_endpoints_enabled") > 0


def test_compute_failed_endpoint_maps_invalid_transition_to_conflict(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        "job-1",
        {"tenant_id": "tenant-worker", "request_id": "req-1", "user_id": "user-1"},
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1"
    ).token

    def _raise_invalid_transition(**_kwargs):
        raise ValueError("invalid job state transition: finished -> queued")

    monkeypatch.setattr(zpe_router, "submit_compute_failure", _raise_invalid_transition)

    response = client.post(
        "/api/zpe/compute/jobs/job-1/failed",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-worker",
            "lease_id": "lease-1",
            "error_code": "ERR",
            "error_message": "boom",
        },
    )

    assert response.status_code == 409
    payload = response.json()
    assert payload["error"]["code"] == "conflict"
    assert "invalid job state transition" in payload["error"]["message"]


def test_compute_failed_endpoint_rejects_tenant_boundary_violation(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        "job-tenant-check",
        {"tenant_id": "tenant-worker", "request_id": "req-1", "user_id": "user-1"},
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1"
    ).token

    response = client.post(
        "/api/zpe/compute/jobs/job-tenant-check/failed",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-other",
            "lease_id": "lease-1",
            "error_code": "ERR",
            "error_message": "boom",
        },
    )

    assert response.status_code == 403
    assert response.json()["error"]["message"] == "tenant boundary violation"


def _result_execution_event(
    *,
    tenant_id: str,
    workspace_id: str,
    job_id: str,
) -> dict[str, object]:
    return {
        "event_id": f"evt-result-{job_id}",
        "tenant_id": tenant_id,
        "workspace_id": workspace_id,
        "job_id": job_id,
        "submission_id": "submission-1",
        "execution_id": "execution-1",
        "occurred_at": "2026-01-01T00:00:00+00:00",
        "trace_id": "trace-1",
        "state": "completed",
        "result_ref": {"output_uri": "zpe://jobs/job-1/result"},
    }


def _failed_execution_event(
    *,
    tenant_id: str,
    workspace_id: str,
    job_id: str,
    error_code: str,
    error_message: str,
) -> dict[str, object]:
    return {
        "event_id": f"evt-failed-{job_id}",
        "tenant_id": tenant_id,
        "workspace_id": workspace_id,
        "job_id": job_id,
        "submission_id": "submission-1",
        "execution_id": "execution-1",
        "occurred_at": "2026-01-01T00:00:01+00:00",
        "trace_id": "trace-2",
        "state": "failed",
        "error": {
            "code": error_code,
            "message": error_message,
            "retryable": True,
        },
    }


def test_compute_result_endpoint_rejects_worker_token_tenant_scope_violation(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    job_id = "job-result-tenant-scope"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": "tenant-worker",
            "workspace_id": "workspace-1",
            "request_id": "req-1",
            "user_id": "user-1",
        },
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1",
        tenant_id="tenant-other",
    ).token

    response = client.post(
        f"/api/zpe/compute/jobs/{job_id}/result",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-worker",
            "lease_id": "lease-1",
            "result": {"zpe_ev": 0.123},
            "summary_text": "ok",
            "freqs_csv": "frequency_cm^-1,intensity\n100,1.0",
        },
    )

    assert response.status_code == 403
    assert response.json()["error"]["message"] == "worker token tenant scope violation"


def test_compute_failed_endpoint_rejects_worker_token_workspace_scope_violation(
    monkeypatch,
):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    job_id = "job-failed-workspace-scope"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": "tenant-worker",
            "workspace_id": "workspace-1",
            "request_id": "req-1",
            "user_id": "user-1",
        },
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1",
        tenant_id="tenant-worker",
        workspace_id="workspace-other",
    ).token

    response = client.post(
        f"/api/zpe/compute/jobs/{job_id}/failed",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-worker",
            "lease_id": "lease-1",
            "error_code": "ERR",
            "error_message": "boom",
        },
    )

    assert response.status_code == 403
    assert (
        response.json()["error"]["message"] == "worker token workspace scope violation"
    )


def test_compute_result_endpoint_rejects_execution_event_workspace_spoof(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    job_id = "job-result-workspace-boundary"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": "tenant-worker",
            "workspace_id": "workspace-1",
            "request_id": "req-1",
            "user_id": "user-1",
        },
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1"
    ).token

    response = client.post(
        f"/api/zpe/compute/jobs/{job_id}/result",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-worker",
            "lease_id": "lease-1",
            "result": {"zpe_ev": 0.123},
            "summary_text": "ok",
            "freqs_csv": "frequency_cm^-1,intensity\n100,1.0",
            "execution_event": _result_execution_event(
                tenant_id="tenant-worker",
                workspace_id="workspace-spoofed",
                job_id=job_id,
            ),
        },
    )

    assert response.status_code == 403
    assert response.json()["error"]["message"] == "workspace boundary violation"


def test_compute_failed_endpoint_rejects_execution_event_workspace_mismatch(monkeypatch):
    fake = _patch_redis(monkeypatch)
    client = TestClient(main.app)
    job_id = "job-failed-workspace-boundary"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": "tenant-worker",
            "workspace_id": "workspace-1",
            "request_id": "req-1",
            "user_id": "user-1",
        },
    )
    worker_token = zpe_worker_auth.WorkerTokenStore(redis=fake).create_token(
        "compute-1"
    ).token

    response = client.post(
        f"/api/zpe/compute/jobs/{job_id}/failed",
        headers={"Authorization": f"Bearer {worker_token}"},
        json={
            "tenant_id": "tenant-worker",
            "lease_id": "lease-1",
            "error_code": "ERR",
            "error_message": "boom",
            "execution_event": _failed_execution_event(
                tenant_id="tenant-worker",
                workspace_id="workspace-other",
                job_id=job_id,
                error_code="ERR",
                error_message="boom",
            ),
        },
    )

    assert response.status_code == 403
    assert response.json()["error"]["message"] == "workspace boundary violation"


def test_compute_failed_endpoint_retries_after_audit_outage(monkeypatch):
    fake = _patch_redis(monkeypatch)
    job_id = "job-audit-retry"
    lease_id = "lease-1"
    worker_id = "compute-1"
    fake.hset(
        f"zpe:status:{job_id}",
        mapping={
            "status": "started",
            "detail": "",
            "updated_at": "2026-01-01T00:00:00+00:00",
        },
    )
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={
            "worker_id": worker_id,
            "lease_id": lease_id,
        },
    )
    monkeypatch.setattr(
        zpe_compute_results, "_dispatch_runtime_state_transition", lambda **_kwargs: None
    )

    first = zpe_compute_results.submit_failure(
        job_id=job_id,
        worker_id=worker_id,
        lease_id=lease_id,
        error_code="ERR_TEMP",
        error_message="temporary outage",
    )
    assert first.requeued is True
    assert first.retry_count == 1
    assert fake.get(f"zpe:retry_count:{job_id}") == b"1"
    assert fake.hget(f"zpe:status:{job_id}", "status") == b"queued"
    assert fake.hgetall(f"zpe:lease:{job_id}") == {}

    second = zpe_compute_results.submit_failure(
        job_id=job_id,
        worker_id=worker_id,
        lease_id=lease_id,
        error_code="ERR_TEMP",
        error_message="temporary outage",
    )
    assert second.requeued is True
    assert second.retry_count == 1
    assert fake.get(f"zpe:retry_count:{job_id}") == b"1"
    assert fake.hget(f"zpe:status:{job_id}", "status") == b"queued"


def test_compute_failed_replay_returns_recorded_outcome_when_lease_already_consumed(
    monkeypatch,
):
    fake = _patch_redis(monkeypatch)
    job_id = "job-overlap-replay"
    lease_id = "lease-overlap"
    worker_id = "compute-1"
    submit_key = f"zpe:failure_submit:{job_id}:{lease_id}"

    fake.hset(
        f"zpe:status:{job_id}",
        mapping={
            "status": "started",
            "detail": "",
            "updated_at": "2026-01-01T00:00:00+00:00",
        },
    )
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": worker_id, "lease_id": lease_id},
    )
    monkeypatch.setattr(
        zpe_compute_results, "_dispatch_runtime_state_transition", lambda **_kwargs: None
    )

    original_get_lease = zpe_compute_results._get_lease
    injected = {"done": False}

    def _get_lease_with_overlap(redis_client, target_job_id):
        if not injected["done"] and target_job_id == job_id:
            injected["done"] = True
            redis_client.hset(
                submit_key,
                mapping={
                    "worker_id": worker_id,
                    "error_code": "ERR_TEMP",
                    "error_message": "temporary outage",
                    "traceback": "",
                    "retry_count": "1",
                    "requeued": "1",
                },
            )
            redis_client.delete(f"zpe:lease:{job_id}")
            return {}
        return original_get_lease(redis_client, target_job_id)

    monkeypatch.setattr(zpe_compute_results, "_get_lease", _get_lease_with_overlap)

    outcome = zpe_compute_results.submit_failure(
        job_id=job_id,
        worker_id=worker_id,
        lease_id=lease_id,
        error_code="ERR_TEMP",
        error_message="temporary outage",
    )

    assert outcome.requeued is True
    assert outcome.retry_count == 1
    assert injected["done"] is True
    assert fake.get(f"zpe:retry_count:{job_id}") is None
    assert fake.hget(f"zpe:status:{job_id}", "status") == b"started"


def test_compute_failed_replay_rejects_same_key_with_different_payload(monkeypatch):
    fake = _patch_redis(monkeypatch)
    job_id = "job-replay-mismatch"
    lease_id = "lease-mismatch"
    worker_id = "compute-1"
    submit_key = f"zpe:failure_submit:{job_id}:{lease_id}"

    fake.hset(
        submit_key,
        mapping={
            "worker_id": worker_id,
            "error_code": "ERR_TEMP",
            "error_message": "temporary outage",
            "traceback": "",
            "retry_count": "1",
            "requeued": "1",
        },
    )

    with pytest.raises(ValueError, match="different payload"):
        zpe_compute_results.submit_failure(
            job_id=job_id,
            worker_id=worker_id,
            lease_id=lease_id,
            error_code="ERR_TEMP",
            error_message="changed message",
        )


def test_zpe_enqueue_returns_503_when_audit_write_fails(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_backends, "enqueue_zpe_job", lambda *_args, **_kwargs: "job-audit")
    monkeypatch.setattr(
        zpe_router, "write_audit_event", lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("disk full"))
    )

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 503
    assert response.json()["error"]["message"] == "audit sink unavailable"
    owner = zpe_job_owner.JobOwnerStore(redis=fake).get_owner_record("job-audit")
    assert owner is not None
    assert owner.user_id == user_id
    assert owner.tenant_id == headers["X-Tenant-Id"]


def test_zpe_retry_after_transient_audit_503_does_not_duplicate_enqueue(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")
    enqueue_calls: list[dict[str, object]] = []

    def _capture_enqueue(payload, *, queue_name=None):
        enqueue_calls.append({"queue_name": queue_name, "payload": payload})
        return "job-audit-retry"

    audit_attempts = {"count": 0}

    def _flaky_audit(**_kwargs):
        audit_attempts["count"] += 1
        if audit_attempts["count"] == 1:
            raise RuntimeError("disk full")
        return False

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_backends, "enqueue_zpe_job", _capture_enqueue)
    monkeypatch.setattr(zpe_router, "write_audit_event", _flaky_audit)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    headers["X-Request-Id"] = "retry-audit-503-request-id"
    payload = {
        "content": QE_INPUT,
        "mobile_indices": [0],
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "continue",
    }

    first = client.post("/api/zpe/jobs", json=payload, headers=headers)
    assert first.status_code == 503
    assert first.json()["error"]["message"] == "audit sink unavailable"

    second = client.post("/api/zpe/jobs", json=payload, headers=headers)
    assert second.status_code == 200
    assert second.json()["id"] == "job-audit-retry"
    assert len(enqueue_calls) == 1


def test_zpe_submit_idempotency_prevents_duplicate_enqueue_under_overlap(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(compute_mode="mock", result_store="redis")
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    tenant_id = headers["X-Tenant-Id"]
    request_payload = {
        "content": QE_INPUT,
        "mobile_indices": [0],
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "continue",
    }
    request_model = ZPEJobRequest(**request_payload)
    request_id = "overlap-request-id-1"
    enqueue_calls = {"count": 0}
    enqueue_lock = threading.Lock()
    barrier = threading.Barrier(2)

    def _slow_enqueue(_payload, *, queue_name=None):
        with enqueue_lock:
            enqueue_calls["count"] += 1
        # Keep claim pending long enough so the overlapping caller takes loser path.
        time.sleep(0.2)
        return "job-overlap-idempotent"

    monkeypatch.setattr(zpe_backends, "enqueue_zpe_job", _slow_enqueue)

    def _submit_once():
        barrier.wait(timeout=1.0)
        return zpe_backends.submit_job(
            request_model,
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [executor.submit(_submit_once) for _ in range(2)]
        outcomes = [future.result(timeout=3.0) for future in futures]

    assert enqueue_calls["count"] == 1
    assert outcomes[0].job_id == "job-overlap-idempotent"
    assert outcomes[1].job_id == "job-overlap-idempotent"
    assert sorted([outcomes[0].idempotent, outcomes[1].idempotent]) == [False, True]


def test_zpe_submit_idempotency_rejects_intent_change_when_target_changes(monkeypatch):
    fake = _patch_redis(monkeypatch)
    settings = ZPESettings(compute_mode="mock", result_store="redis")
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_router, "get_zpe_settings", lambda: settings)
    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(
        client, monkeypatch, fake, queue_name="queue-a"
    )
    headers["X-Request-Id"] = "idempotency-intent-change-1"

    enroll = client.post(
        "/api/zpe/compute/enroll-tokens",
        json={"ttl_seconds": 60},
        headers=headers,
    )
    assert enroll.status_code == 200
    enroll_token = enroll.json()["token"]
    register = client.post(
        "/api/zpe/compute/servers",
        json={
            "token": enroll_token,
            "name": "server-2",
            "queue_name": "queue-b",
            "activate_target": False,
        },
    )
    assert register.status_code == 200

    request_payload = {
        "content": QE_INPUT,
        "mobile_indices": [0],
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "continue",
    }
    first = client.post("/api/zpe/jobs", json=request_payload, headers=headers)
    assert first.status_code == 200

    targets = client.get("/api/zpe/targets", headers=headers)
    assert targets.status_code == 200
    queue_b_target_id = next(
        target["id"]
        for target in targets.json()["targets"]
        if target["queue_name"] == "queue-b"
    )
    switched = client.put(
        f"/api/zpe/targets/{queue_b_target_id}/active",
        headers=headers,
    )
    assert switched.status_code == 200

    second = client.post("/api/zpe/jobs", json=request_payload, headers=headers)
    assert second.status_code == 409
    assert "different submission payload" in second.json()["error"]["message"]


def test_submit_idempotency_claim_uses_legacy_ready_record_during_key_migration(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = zpe_submit_idempotency.SubmitIdempotencyStore(redis=fake)
    payload = {
        "content": QE_INPUT,
        "mobile_indices": [0],
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "continue",
    }
    request_fingerprint = zpe_submit_idempotency.compute_submit_request_fingerprint(payload)
    tenant_id = "tenant-legacy"
    user_id = "user-legacy"
    request_id = "legacy-request-id-1"
    # Simulate a pre-migration key written as tenant_id:user_id:request_id.
    legacy_key = f"zpe:submit:idempotency:{tenant_id}:{user_id}:{request_id}"
    fake.set(
        legacy_key,
        json.dumps(
            {
                "state": "ready",
                "job_id": "job-from-legacy-key",
                "request_fingerprint": request_fingerprint,
                "requested_queue_name": "zpe",
                "resolved_queue_name": "zpe",
            },
            ensure_ascii=True,
            separators=(",", ":"),
        ),
    )

    claim = store.claim_or_get(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        request_fingerprint=request_fingerprint,
    )

    assert claim.state == "ready"
    assert claim.record is not None
    assert claim.record.job_id == "job-from-legacy-key"


def test_submit_idempotency_claim_rejects_payload_change_for_legacy_key(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = zpe_submit_idempotency.SubmitIdempotencyStore(redis=fake)
    tenant_id = "tenant-legacy"
    user_id = "user-legacy"
    request_id = "legacy-request-id-2"
    legacy_key = f"zpe:submit:idempotency:{tenant_id}:{user_id}:{request_id}"
    fake.set(
        legacy_key,
        json.dumps(
            {
                "state": "pending",
                "request_fingerprint": "old-fingerprint",
                "claim_token": "legacy-claim",
            },
            ensure_ascii=True,
            separators=(",", ":"),
        ),
    )

    with pytest.raises(ValueError, match="different submission payload"):
        store.claim_or_get(
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
            request_fingerprint="new-fingerprint",
        )


def test_zpe_owner_enforcement_blocks_non_owner(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    owner_headers, _owner_id = _setup_user_and_target(
        client, monkeypatch, fake, user_id="user-owner"
    )
    other_headers = {
        "X-Dev-User-Id": "user-other",
        "X-Dev-User-Email": "user-other@example.com",
        "X-Tenant-Id": "tenant-user-other",
    }

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=owner_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    status_forbidden = client.get(f"/api/zpe/jobs/{job_id}", headers=other_headers)
    assert status_forbidden.status_code == 403

    result_forbidden = client.get(
        f"/api/zpe/jobs/{job_id}/result", headers=other_headers
    )
    assert result_forbidden.status_code == 403

    file_forbidden = client.get(
        f"/api/zpe/jobs/{job_id}/files",
        params={"kind": "summary"},
        headers=other_headers,
    )
    assert file_forbidden.status_code == 403


def test_zpe_owner_enforcement_blocks_wrong_tenant(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    owner_headers, _owner_id = _setup_user_and_target(
        client, monkeypatch, fake, user_id="tenant-owner"
    )
    wrong_tenant_headers = {
        "X-Dev-User-Id": "tenant-owner",
        "X-Dev-User-Email": "tenant-owner@example.com",
        "X-Tenant-Id": "tenant-other",
    }

    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=owner_headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    forbidden = client.get(f"/api/zpe/jobs/{job_id}", headers=wrong_tenant_headers)
    assert forbidden.status_code == 403


def test_zpe_owner_enforcement_prefers_job_meta_over_legacy_owner_key(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")
    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, user_id = _setup_user_and_target(client, monkeypatch, fake)
    job_id = "job-owner-priority"
    zpe_job_meta.JobMetaStore(redis=fake).set_meta(
        job_id,
        {
            "tenant_id": headers["X-Tenant-Id"],
            "user_id": user_id,
            "request_id": "req-owner-priority",
        },
    )
    fake.set(
        f"zpe:job:owner:{job_id}",
        json.dumps(
            {"user_id": "legacy-other-user", "tenant_id": "tenant-legacy-other"},
            ensure_ascii=True,
            separators=(",", ":"),
        ),
    )
    store.set_status(job_id, "queued")

    allowed = client.get(f"/api/zpe/jobs/{job_id}", headers=headers)
    assert allowed.status_code == 200


def test_zpe_requires_tenant_id(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    headers.pop("X-Tenant-Id")
    response = client.post(
        "/api/zpe/jobs",
        json={
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 400
    assert response.json()["error"]["message"] == "missing tenant_id"


def test_qe_relax_v1_enqueue_run_and_result(monkeypatch):
    fake = _patch_redis(monkeypatch)
    store = RedisResultStore(redis=fake)
    settings = ZPESettings(compute_mode="mock", result_store="redis")

    monkeypatch.setattr(zpe_backends, "get_result_store", lambda: store)
    monkeypatch.setattr(zpe_backends, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    headers, _user_id = _setup_user_and_target(client, monkeypatch, fake)
    response = client.post(
        "/api/zpe/jobs",
        json={
            "calc_type": "qe.relax.v1",
            "content": QE_INPUT,
            "mobile_indices": [0],
            "use_environ": False,
            "input_dir": None,
            "calc_mode": "continue",
        },
        headers=headers,
    )
    assert response.status_code == 200
    job_id = response.json()["id"]

    status = client.get(f"/api/zpe/jobs/{job_id}", headers=headers)
    assert status.status_code == 200
    assert status.json()["status"] == "finished"

    result = client.get(f"/api/zpe/jobs/{job_id}/result", headers=headers)
    assert result.status_code == 200
    payload = result.json()["result"]
    assert payload["calc_type"] == "qe.relax.v1"
    assert isinstance(payload["freqs_cm"], list)
