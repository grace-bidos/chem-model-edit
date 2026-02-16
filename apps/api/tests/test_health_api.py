from __future__ import annotations

from email.message import Message
from io import BytesIO
import subprocess
from typing import Any, Literal
from urllib import error as urlerror

from fastapi.testclient import TestClient

import main


_AIIA_ENV_KEYS = [
    "AIIA_MANAGED_CHECKS_ENABLED",
    "AIIA_MANAGED_BASE_URL",
    "AIIA_MANAGED_HEALTH_PATH",
    "AIIA_MANAGED_READY_PATH",
    "AIIA_MANAGED_TIMEOUT_SECONDS",
    "AIIA_MANAGED_BEARER_TOKEN",
    "AIIA_USER_MANAGED_DEEP_READY_ENABLED",
    "AIIA_USER_MANAGED_DEEP_READY_TIMEOUT_SECONDS",
    "AIIDA_PROFILE",
]


def _clear_aiida_env(monkeypatch: Any) -> None:
    for key in _AIIA_ENV_KEYS:
        monkeypatch.delenv(key, raising=False)


class _DummyResponse:
    def __init__(self, status: int) -> None:
        self.status = status

    def __enter__(self) -> "_DummyResponse":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: Any,
    ) -> Literal[False]:
        _ = (exc_type, exc, traceback)
        return False


def test_health_reports_skipped_when_managed_aiida_checks_disabled(monkeypatch: Any) -> None:
    _clear_aiida_env(monkeypatch)
    client = TestClient(main.app)

    response = client.get("/api/health")

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "ok"
    assert payload["checks"]["managed_aiida_runtime"]["status"] == "skipped"
    assert payload["checks"]["user_managed_deep_readiness"]["status"] == "skipped"


def test_ready_returns_503_when_managed_aiida_is_misconfigured(monkeypatch: Any) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert payload["status"] == "not_ready"
    assert check["status"] == "failed"
    assert check["error"] == "misconfigured"
    assert payload["checks"]["user_managed_deep_readiness"]["status"] == "skipped"


def test_ready_returns_503_when_managed_aiida_base_url_is_malformed(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    monkeypatch.setenv("AIIA_MANAGED_BASE_URL", "http://[::1")
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert payload["status"] == "not_ready"
    assert check["status"] == "failed"
    assert check["error"] == "misconfigured"
    assert "invalid AIIA_MANAGED_BASE_URL" in check["detail"]


def test_health_returns_200_degraded_when_managed_aiida_base_url_is_malformed(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    monkeypatch.setenv("AIIA_MANAGED_BASE_URL", "http://[::1")
    client = TestClient(main.app)

    response = client.get("/api/health")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert payload["status"] == "degraded"
    assert check["status"] == "failed"
    assert check["error"] == "misconfigured"
    assert "invalid AIIA_MANAGED_BASE_URL" in check["detail"]


def test_ready_returns_503_when_managed_aiida_is_unreachable(monkeypatch: Any) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    monkeypatch.setenv("AIIA_MANAGED_BASE_URL", "https://aiida.example")

    def _raise_unreachable(request: Any, timeout: int = 0) -> Any:
        _ = (request, timeout)
        raise urlerror.URLError("connection refused")

    monkeypatch.setattr(
        "services.aiida_runtime_checks.urlrequest.urlopen",
        _raise_unreachable,
    )
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert check["status"] == "failed"
    assert check["error"] == "unreachable"


def test_ready_returns_200_when_managed_aiida_is_healthy(monkeypatch: Any) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    monkeypatch.setenv("AIIA_MANAGED_BASE_URL", "https://aiida.example")
    monkeypatch.setenv("AIIA_MANAGED_READY_PATH", "/api/v1/readyz")

    def _return_healthy(request: Any, timeout: int = 0) -> _DummyResponse:
        _ = timeout
        assert request.full_url == "https://aiida.example/api/v1/readyz"
        return _DummyResponse(status=204)

    monkeypatch.setattr(
        "services.aiida_runtime_checks.urlrequest.urlopen",
        _return_healthy,
    )
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert payload["status"] == "ready"
    assert check["status"] == "ok"
    assert check["http_status"] == 204
    assert payload["checks"]["user_managed_deep_readiness"]["status"] == "skipped"


def test_health_reports_degraded_when_managed_aiida_returns_http_error(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_MANAGED_CHECKS_ENABLED", "true")
    monkeypatch.setenv("AIIA_MANAGED_BASE_URL", "https://aiida.example")

    def _raise_http_error(request: Any, timeout: int = 0) -> Any:
        _ = (request, timeout)
        raise urlerror.HTTPError(
            url="https://aiida.example/health",
            code=503,
            msg="service unavailable",
            hdrs=Message(),
            fp=BytesIO(b""),
        )

    monkeypatch.setattr(
        "services.aiida_runtime_checks.urlrequest.urlopen",
        _raise_http_error,
    )
    client = TestClient(main.app)

    response = client.get("/api/health")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["managed_aiida_runtime"]
    assert payload["status"] == "degraded"
    assert check["status"] == "failed"
    assert check["error"] == "upstream_unhealthy"
    assert check["http_status"] == 503


def test_ready_returns_503_when_user_managed_deep_readiness_tooling_missing(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        if command[0] == "scontrol":
            raise FileNotFoundError("scontrol")
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="ok",
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "not_ready"
    assert check["status"] == "failed"
    assert check["error"] == "tooling_missing"
    assert check["failed_probe"] == "scontrol_ping"
    assert check["probes"]["scontrol_ping"]["error"] == "tooling_missing"


def test_health_reports_degraded_when_user_managed_deep_readiness_command_fails(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")
    commands: list[tuple[str, ...]] = []

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        commands.append(tuple(command))
        if command[0] == "sinfo":
            return subprocess.CompletedProcess(
                args=command,
                returncode=1,
                stdout="",
                stderr="slurmctld: down",
            )
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="ok",
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/health")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "degraded"
    assert check["status"] == "failed"
    assert check["error"] == "command_failed"
    assert check["failed_probe"] == "sinfo"
    assert check["probes"]["sinfo"]["exit_code"] == 1
    assert commands == [("scontrol", "ping"), ("sinfo",)]


def test_health_reports_degraded_when_user_managed_deep_readiness_command_launch_fails(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        if command[0] == "sinfo":
            raise PermissionError("permission denied")
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout="ok",
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/health")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "degraded"
    assert check["status"] == "failed"
    assert check["error"] == "tooling_failure"
    assert check["failed_probe"] == "sinfo"
    assert check["probes"]["sinfo"]["error"] == "tooling_failure"


def test_ready_returns_200_when_user_managed_deep_readiness_passes(monkeypatch: Any) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")
    monkeypatch.setenv("AIIDA_PROFILE", "default")

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        stdout = "ok"
        if command[:3] == ["verdi", "profile", "list"]:
            stdout = "* default\n  other"
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=stdout,
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "ready"
    assert check["status"] == "ok"
    assert set(check["probes"]) == {
        "scontrol_ping",
        "sinfo",
        "verdi_status",
        "verdi_profile",
    }


def test_ready_returns_503_when_aiida_profile_is_missing_from_verdi_profile_list(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")
    monkeypatch.setenv("AIIDA_PROFILE", "target-profile")

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        stdout = "ok"
        if command[:3] == ["verdi", "profile", "list"]:
            stdout = "* default\n  staging"
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=stdout,
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "not_ready"
    assert check["status"] == "failed"
    assert check["error"] == "misconfigured"
    assert check["failed_probe"] == "verdi_profile"


def test_ready_returns_503_when_aiida_profile_matches_only_as_substring(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")
    monkeypatch.setenv("AIIDA_PROFILE", "fault")

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        stdout = "ok"
        if command[:3] == ["verdi", "profile", "list"]:
            stdout = "* default\n  staging"
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=stdout,
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 503
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "not_ready"
    assert check["status"] == "failed"
    assert check["error"] == "misconfigured"
    assert check["failed_probe"] == "verdi_profile"


def test_ready_returns_200_when_aiida_profile_exists_beyond_truncation_boundary(
    monkeypatch: Any,
) -> None:
    _clear_aiida_env(monkeypatch)
    monkeypatch.setenv("AIIA_USER_MANAGED_DEEP_READY_ENABLED", "true")
    monkeypatch.setenv("AIIDA_PROFILE", "target-profile")

    long_profile_list = "\n".join([f"  profile-{idx:03d}" for idx in range(45)])
    long_profile_list = f"{long_profile_list}\n* target-profile"

    def _run_probe(command: list[str], **kwargs: Any) -> Any:
        _ = kwargs
        stdout = "ok"
        if command[:3] == ["verdi", "profile", "list"]:
            stdout = long_profile_list
        return subprocess.CompletedProcess(
            args=command,
            returncode=0,
            stdout=stdout,
            stderr="",
        )

    monkeypatch.setattr("services.aiida_runtime_checks.subprocess.run", _run_probe)
    client = TestClient(main.app)

    response = client.get("/api/ready")

    assert response.status_code == 200
    payload = response.json()
    check = payload["checks"]["user_managed_deep_readiness"]
    assert payload["status"] == "ready"
    assert check["status"] == "ok"
