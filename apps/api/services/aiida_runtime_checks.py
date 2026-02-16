from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import subprocess
import time
from typing import Any, Literal
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

from services.zpe.settings import ZPESettings
from services.zpe.slurm_policy import (
    SlurmPolicyConfigError,
    parse_slurm_adapter_mode,
    parse_slurm_adapter_rollback_guard,
    resolve_effective_slurm_adapter_mode,
    validate_slurm_policy_file,
)


ManagedAiidaStatus = Literal["ok", "failed", "skipped"]

_ENABLED_KEY = "AIIA_MANAGED_CHECKS_ENABLED"
_BASE_URL_KEY = "AIIA_MANAGED_BASE_URL"
_HEALTH_PATH_KEY = "AIIA_MANAGED_HEALTH_PATH"
_READY_PATH_KEY = "AIIA_MANAGED_READY_PATH"
_TIMEOUT_SECONDS_KEY = "AIIA_MANAGED_TIMEOUT_SECONDS"
_TOKEN_KEY = "AIIA_MANAGED_BEARER_TOKEN"
_DEEP_READY_ENABLED_KEY = "AIIA_USER_MANAGED_DEEP_READY_ENABLED"
_DEEP_READY_TIMEOUT_SECONDS_KEY = "AIIA_USER_MANAGED_DEEP_READY_TIMEOUT_SECONDS"
_AIIDA_PROFILE_KEY = "AIIDA_PROFILE"
_SLURM_POLICY_PATH_KEY = "ZPE_SLURM_POLICY_PATH"
_MAX_OUTPUT_LENGTH = 240


def _load_slurm_adapter_config() -> tuple[str | None, str | None, str]:
    """Load adapter config via ZPESettings so .env-backed values are respected."""
    settings = ZPESettings()
    return (
        settings.slurm_adapter,
        settings.slurm_adapter_rollback_guard,
        (settings.slurm_policy_path or "").strip(),
    )


def _truthy(value: str | None) -> bool:
    if value is None:
        return False
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _timeout_seconds() -> int:
    raw = os.getenv(_TIMEOUT_SECONDS_KEY, "3").strip()
    try:
        parsed = int(raw)
    except ValueError:
        return 3
    return parsed if parsed > 0 else 3


def _deep_timeout_seconds() -> int:
    raw = os.getenv(_DEEP_READY_TIMEOUT_SECONDS_KEY, "5").strip()
    try:
        parsed = int(raw)
    except ValueError:
        return 5
    return parsed if parsed > 0 else 5


def _normalize_output(value: str | bytes | None) -> str | None:
    if value is None:
        return None
    if isinstance(value, bytes):
        text = value.decode("utf-8", errors="replace")
    else:
        text = value
    normalized = text.strip()
    if not normalized:
        return None
    if len(normalized) > _MAX_OUTPUT_LENGTH:
        return f"{normalized[:_MAX_OUTPUT_LENGTH]}..."
    return normalized


@dataclass(frozen=True)
class ManagedAiidaCheck:
    status: ManagedAiidaStatus
    detail: str
    url: str | None = None
    error: str | None = None
    http_status: int | None = None
    latency_ms: int | None = None

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "status": self.status,
            "detail": self.detail,
        }
        if self.url is not None:
            payload["url"] = self.url
        if self.error is not None:
            payload["error"] = self.error
        if self.http_status is not None:
            payload["http_status"] = self.http_status
        if self.latency_ms is not None:
            payload["latency_ms"] = self.latency_ms
        return payload


@dataclass(frozen=True)
class CommandProbe:
    status: ManagedAiidaStatus
    detail: str
    command: str
    error: str | None = None
    exit_code: int | None = None
    latency_ms: int | None = None
    stdout: str | None = None
    stderr: str | None = None
    raw_stdout: str | None = None
    raw_stderr: str | None = None

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "status": self.status,
            "detail": self.detail,
            "command": self.command,
        }
        if self.error is not None:
            payload["error"] = self.error
        if self.exit_code is not None:
            payload["exit_code"] = self.exit_code
        if self.latency_ms is not None:
            payload["latency_ms"] = self.latency_ms
        if self.stdout is not None:
            payload["stdout"] = self.stdout
        if self.stderr is not None:
            payload["stderr"] = self.stderr
        return payload


@dataclass(frozen=True)
class UserManagedDeepReadinessCheck:
    status: ManagedAiidaStatus
    detail: str
    probes: dict[str, CommandProbe]
    error: str | None = None
    failed_probe: str | None = None

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "status": self.status,
            "detail": self.detail,
            "probes": {
                probe_name: probe_result.as_dict()
                for probe_name, probe_result in self.probes.items()
            },
        }
        if self.error is not None:
            payload["error"] = self.error
        if self.failed_probe is not None:
            payload["failed_probe"] = self.failed_probe
        return payload


@dataclass(frozen=True)
class SlurmRealAdapterPreconditionCheck:
    status: ManagedAiidaStatus
    detail: str
    adapter_configured: str
    adapter_effective: str
    rollback_guard: str
    probes: dict[str, CommandProbe]
    error: str | None = None
    failed_probe: str | None = None

    def as_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "status": self.status,
            "detail": self.detail,
            "adapter_configured": self.adapter_configured,
            "adapter_effective": self.adapter_effective,
            "rollback_guard": self.rollback_guard,
            "probes": {
                probe_name: probe_result.as_dict()
                for probe_name, probe_result in self.probes.items()
            },
        }
        if self.error is not None:
            payload["error"] = self.error
        if self.failed_probe is not None:
            payload["failed_probe"] = self.failed_probe
        return payload


def _run_command_probe(*, command: list[str]) -> CommandProbe:
    started = time.perf_counter()
    command_text = " ".join(command)
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=_deep_timeout_seconds(),
        )
    except FileNotFoundError:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return CommandProbe(
            status="failed",
            detail=f"required command missing: {command[0]}",
            command=command_text,
            error="tooling_missing",
            latency_ms=latency_ms,
        )
    except OSError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return CommandProbe(
            status="failed",
            detail=f"command launch failed: {command_text}",
            command=command_text,
            error="tooling_failure",
            latency_ms=latency_ms,
            stderr=_normalize_output(str(exc)),
            raw_stderr=str(exc),
        )
    except subprocess.TimeoutExpired as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return CommandProbe(
            status="failed",
            detail=f"command timed out after {_deep_timeout_seconds()}s: {command_text}",
            command=command_text,
            error="timeout",
            latency_ms=latency_ms,
            stdout=_normalize_output(exc.stdout),
            stderr=_normalize_output(exc.stderr),
            raw_stdout=_normalize_output(exc.stdout),
            raw_stderr=_normalize_output(exc.stderr),
        )

    latency_ms = int((time.perf_counter() - started) * 1000)
    raw_stdout = completed.stdout
    raw_stderr = completed.stderr
    stdout = _normalize_output(completed.stdout)
    stderr = _normalize_output(completed.stderr)
    if completed.returncode != 0:
        return CommandProbe(
            status="failed",
            detail=f"command failed with exit {completed.returncode}: {command_text}",
            command=command_text,
            error="command_failed",
            exit_code=int(completed.returncode),
            latency_ms=latency_ms,
            stdout=stdout,
            stderr=stderr,
            raw_stdout=raw_stdout,
            raw_stderr=raw_stderr,
        )

    return CommandProbe(
        status="ok",
        detail="command succeeded",
        command=command_text,
        latency_ms=latency_ms,
        stdout=stdout,
        raw_stdout=raw_stdout,
        raw_stderr=raw_stderr,
    )


def _extract_profile_names(output: str | None) -> set[str]:
    names: set[str] = set()
    if output is None:
        return names

    for raw_line in output.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("*"):
            line = line[1:].strip()
        if not line:
            continue
        names.add(line.split()[0])
    return names


def probe_user_managed_deep_readiness(
    *, check_kind: Literal["health", "ready"]
) -> UserManagedDeepReadinessCheck:
    _ = check_kind
    enabled = _truthy(os.getenv(_DEEP_READY_ENABLED_KEY))
    if not enabled:
        return UserManagedDeepReadinessCheck(
            status="skipped",
            detail="user-managed deep readiness checks disabled",
            probes={},
        )

    probe_specs: list[tuple[str, list[str]]] = [
        ("scontrol_ping", ["scontrol", "ping"]),
        ("sinfo", ["sinfo"]),
        ("verdi_status", ["verdi", "status"]),
        ("verdi_profile", ["verdi", "profile", "list"]),
    ]
    probes: dict[str, CommandProbe] = {}
    for probe_name, command in probe_specs:
        probe = _run_command_probe(command=command)
        probes[probe_name] = probe
        if probe.status == "failed":
            return UserManagedDeepReadinessCheck(
                status="failed",
                detail="user-managed deep readiness probes failed",
                probes=probes,
                error=probe.error,
                failed_probe=probe_name,
            )

    configured_profile = (os.getenv(_AIIDA_PROFILE_KEY) or "").strip()
    verdi_profile_probe = probes.get("verdi_profile")
    if (
        verdi_profile_probe is not None
        and configured_profile
        and verdi_profile_probe.status == "ok"
        and configured_profile not in _extract_profile_names(verdi_profile_probe.raw_stdout)
    ):
        probes["verdi_profile"] = CommandProbe(
            status="failed",
            detail=f"{_AIIDA_PROFILE_KEY} is not present in verdi profile list",
            command=verdi_profile_probe.command,
            error="misconfigured",
            latency_ms=verdi_profile_probe.latency_ms,
            stdout=verdi_profile_probe.stdout,
            stderr=verdi_profile_probe.stderr,
            raw_stdout=verdi_profile_probe.raw_stdout,
            raw_stderr=verdi_profile_probe.raw_stderr,
        )

        return UserManagedDeepReadinessCheck(
            status="failed",
            detail="user-managed deep readiness probes failed",
            probes=probes,
            error="misconfigured",
            failed_probe="verdi_profile",
        )

    return UserManagedDeepReadinessCheck(
        status="ok",
        detail="user-managed deep readiness probes passed",
        probes=probes,
    )


def _policy_probe(policy_path_value: str) -> CommandProbe:
    if not policy_path_value:
        return CommandProbe(
            status="failed",
            detail=f"{_SLURM_POLICY_PATH_KEY} must be set for real-policy",
            command="validate slurm policy",
            error="misconfigured",
        )
    try:
        validate_slurm_policy_file(policy_path=Path(policy_path_value))
    except SlurmPolicyConfigError as exc:
        return CommandProbe(
            status="failed",
            detail=f"invalid slurm policy: {exc}",
            command="validate slurm policy",
            error="misconfigured",
        )

    return CommandProbe(
        status="ok",
        detail="slurm policy file is valid",
        command="validate slurm policy",
    )


def probe_slurm_real_adapter_preconditions(
    *, check_kind: Literal["health", "ready"]
) -> SlurmRealAdapterPreconditionCheck:
    _ = check_kind
    adapter_value, guard_value, policy_path_value = _load_slurm_adapter_config()
    try:
        configured = parse_slurm_adapter_mode(adapter_value)
        guard = parse_slurm_adapter_rollback_guard(guard_value)
        effective = resolve_effective_slurm_adapter_mode(configured, rollback_guard=guard)
    except SlurmPolicyConfigError as exc:
        return SlurmRealAdapterPreconditionCheck(
            status="failed",
            detail=f"real-adapter configuration invalid: {exc}",
            adapter_configured=(adapter_value or "stub-policy"),
            adapter_effective="unknown",
            rollback_guard=(guard_value or "allow"),
            probes={},
            error="misconfigured",
        )

    if configured != "real-policy":
        return SlurmRealAdapterPreconditionCheck(
            status="skipped",
            detail="slurm adapter is not configured as real-policy",
            adapter_configured=configured,
            adapter_effective=effective,
            rollback_guard=guard,
            probes={},
        )

    if effective != "real-policy":
        return SlurmRealAdapterPreconditionCheck(
            status="failed",
            detail="real-policy configured but rollback guard forces fallback adapter",
            adapter_configured=configured,
            adapter_effective=effective,
            rollback_guard=guard,
            probes={},
            error="rollback_guard_active",
        )

    probes: dict[str, CommandProbe] = {}
    policy_probe = _policy_probe(policy_path_value)
    probes["policy_file"] = policy_probe
    if policy_probe.status == "failed":
        return SlurmRealAdapterPreconditionCheck(
            status="failed",
            detail="real-policy precondition checks failed",
            adapter_configured=configured,
            adapter_effective=effective,
            rollback_guard=guard,
            probes=probes,
            error=policy_probe.error,
            failed_probe="policy_file",
        )

    for probe_name, command in (
        ("scontrol_ping", ["scontrol", "ping"]),
        ("sinfo", ["sinfo"]),
    ):
        probe = _run_command_probe(command=command)
        probes[probe_name] = probe
        if probe.status == "failed":
            return SlurmRealAdapterPreconditionCheck(
                status="failed",
                detail="real-policy precondition checks failed",
                adapter_configured=configured,
                adapter_effective=effective,
                rollback_guard=guard,
                probes=probes,
                error=probe.error,
                failed_probe=probe_name,
            )

    return SlurmRealAdapterPreconditionCheck(
        status="ok",
        detail="real-policy precondition checks passed",
        adapter_configured=configured,
        adapter_effective=effective,
        rollback_guard=guard,
        probes=probes,
    )


def probe_managed_aiida_runtime(*, check_kind: Literal["health", "ready"]) -> ManagedAiidaCheck:
    """Probe managed AiiDA runtime according to the selected check kind."""

    enabled = _truthy(os.getenv(_ENABLED_KEY))
    if not enabled:
        return ManagedAiidaCheck(
            status="skipped",
            detail="managed AiiDA checks disabled",
        )

    base_url = (os.getenv(_BASE_URL_KEY) or "").strip().rstrip("/")
    if not base_url:
        return ManagedAiidaCheck(
            status="failed",
            detail=f"{_BASE_URL_KEY} is required when {_ENABLED_KEY}=true",
            error="misconfigured",
        )

    if check_kind == "health":
        path = os.getenv(_HEALTH_PATH_KEY, "/health")
    else:
        path = os.getenv(_READY_PATH_KEY, "/ready")

    try:
        parsed_base = urlparse.urlsplit(base_url)
    except ValueError as exc:
        return ManagedAiidaCheck(
            status="failed",
            detail=f"invalid {_BASE_URL_KEY}: {exc}",
            error="misconfigured",
        )
    if not parsed_base.scheme or not parsed_base.netloc:
        return ManagedAiidaCheck(
            status="failed",
            detail=f"invalid {_BASE_URL_KEY}: scheme and host are required",
            error="misconfigured",
        )

    target_url = f"{base_url}/{path.lstrip('/')}"

    headers = {"Accept": "application/json"}
    token = (os.getenv(_TOKEN_KEY) or "").strip()
    if token:
        headers["Authorization"] = f"Bearer {token}"
    try:
        request = urlrequest.Request(target_url, headers=headers, method="GET")
    except ValueError as exc:
        return ManagedAiidaCheck(
            status="failed",
            detail=f"invalid managed AiiDA probe URL: {exc}",
            url=target_url,
            error="misconfigured",
        )

    timeout = _timeout_seconds()
    started = time.perf_counter()
    try:
        with urlrequest.urlopen(request, timeout=timeout) as response:
            http_status = int(response.status)
            latency_ms = int((time.perf_counter() - started) * 1000)
            return ManagedAiidaCheck(
                status="ok",
                detail="managed AiiDA runtime reachable",
                url=target_url,
                http_status=http_status,
                latency_ms=latency_ms,
            )
    except urlerror.HTTPError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return ManagedAiidaCheck(
            status="failed",
            detail=f"managed AiiDA returned HTTP {exc.code}",
            url=target_url,
            error="upstream_unhealthy",
            http_status=int(exc.code),
            latency_ms=latency_ms,
        )
    except ValueError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return ManagedAiidaCheck(
            status="failed",
            detail=f"invalid managed AiiDA probe URL: {exc}",
            url=target_url,
            error="misconfigured",
            latency_ms=latency_ms,
        )
    except urlerror.URLError as exc:
        latency_ms = int((time.perf_counter() - started) * 1000)
        return ManagedAiidaCheck(
            status="failed",
            detail=f"managed AiiDA endpoint unreachable: {exc.reason}",
            url=target_url,
            error="unreachable",
            latency_ms=latency_ms,
        )
