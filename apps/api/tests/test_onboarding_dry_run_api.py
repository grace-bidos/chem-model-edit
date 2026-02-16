from __future__ import annotations

from fastapi.testclient import TestClient

import main
from app.api import app as api_app
from app import deps as app_deps
from services.zpe.settings import ZPESettings


def _base_policy(*, fallback_mode: str = "route-default") -> dict:
    fallback_policy: dict[str, str] = {
        "mode": fallback_mode,
        "default_queue": "standard",
    }
    if fallback_mode == "deny":
        fallback_policy = {"mode": "deny"}

    return {
        "cluster_name": "chem-byo-cluster-a",
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
        "fallback_policy": fallback_policy,
        "registration_requirements": {
            "required_node_labels": ["owner", "region"],
            "required_health_checks": ["sinfo", "sbatch"],
            "deny_if_drain_or_down": True,
        },
    }


def _base_registration() -> dict:
    return {
        "node_id": "node-01",
        "cluster_name": "chem-byo-cluster-a",
        "enabled_queues": ["standard"],
        "node_labels": {
            "owner": "chem-team",
            "region": "us-central1",
        },
        "health_checks": {
            "sinfo": True,
            "sbatch": True,
        },
        "slurm_state": "idle",
    }


def test_onboarding_dry_run_requires_admin_token(monkeypatch) -> None:
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(app_deps, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/admin/onboarding/dry-run",
        json={
            "path_mode": "path-a-existing-slurm",
            "policy": _base_policy(),
        },
    )

    assert response.status_code == 401


def test_onboarding_route_available_from_app_api_import(monkeypatch) -> None:
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(app_deps, "get_zpe_settings", lambda: settings)

    client = TestClient(api_app)
    response = client.post(
        "/api/zpe/admin/onboarding/dry-run",
        headers={"Authorization": "Bearer secret"},
        json={
            "path_mode": "path-a-existing-slurm",
            "policy": _base_policy(),
        },
    )

    assert response.status_code == 200


def test_onboarding_dry_run_returns_valid_report_and_queue_decision(monkeypatch) -> None:
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(app_deps, "get_zpe_settings", lambda: settings)

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/admin/onboarding/dry-run",
        headers={"Authorization": "Bearer secret"},
        json={
            "path_mode": "path-a-existing-slurm",
            "policy": _base_policy(),
            "registration": _base_registration(),
            "requested_queue": "experimental",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "valid"
    assert payload["errors"] == []
    assert payload["normalized_manifest"]["schema_version"] == "byo_onboarding_manifest/v1"
    assert payload["normalized_manifest"]["path_mode"] == "path-a-existing-slurm"
    queue_decision = payload["report"]["queue_decision"]
    assert queue_decision["requested_queue"] == "experimental"
    assert queue_decision["resolved_queue"] == "standard"
    assert queue_decision["used_fallback"] is True


def test_onboarding_dry_run_returns_structured_errors(monkeypatch) -> None:
    settings = ZPESettings(admin_token="secret")
    monkeypatch.setattr(app_deps, "get_zpe_settings", lambda: settings)

    registration = _base_registration()
    registration["cluster_name"] = "other-cluster"
    registration["health_checks"]["sbatch"] = False

    client = TestClient(main.app)
    response = client.post(
        "/api/zpe/admin/onboarding/dry-run",
        headers={"Authorization": "Bearer secret"},
        json={
            "path_mode": "path-b-bootstrap-slurm-aiida",
            "policy": _base_policy(fallback_mode="deny"),
            "registration": registration,
            "requested_queue": "unknown",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "invalid"
    assert payload["report"]["error_count"] >= 3
    assert payload["report"]["registration_valid"] is False
    assert payload["report"]["queue_resolution_valid"] is False

    sections = {issue["section"] for issue in payload["errors"]}
    assert "registration" in sections
    assert "queue_resolution" in sections

    issue_paths = {issue["field_path"] for issue in payload["errors"]}
    assert "cluster_name" in issue_paths
    assert "health_checks.sbatch" in issue_paths
    assert "requested_queue" in issue_paths

    assert all(issue["action"] for issue in payload["errors"])
