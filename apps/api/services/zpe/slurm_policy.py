from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
import subprocess
from typing import Any, Literal, Mapping, cast


SlurmAdapterMode = Literal["passthrough", "stub-policy", "real-policy"]
SlurmAdapterRollbackGuard = Literal["allow", "force-stub-policy", "force-passthrough"]


@dataclass(frozen=True)
class SlurmQueueMapping:
    queue: str
    partition: str
    account: str
    qos: str
    max_walltime_minutes: int


@dataclass(frozen=True)
class SlurmQueueResolution:
    requested_queue: str
    resolved_queue: str
    used_fallback: bool
    mapping: SlurmQueueMapping


SLURM_ADAPTER_CONTRACT_VERSION = "slurm-adapter-boundary/v1"


@dataclass(frozen=True)
class SlurmAdapterBoundaryStub:
    adapter: SlurmAdapterMode
    contract_version: str
    requested_queue: str
    resolved_queue: str
    used_fallback: bool
    mapping: SlurmQueueMapping | None


class SlurmPolicyError(RuntimeError):
    """Base exception for runtime Slurm policy resolution errors."""


class SlurmPolicyDeniedError(SlurmPolicyError):
    """Raised when policy denies the requested queue."""


class SlurmPolicyConfigError(SlurmPolicyError):
    """Raised when policy content is missing or invalid for runtime resolution."""


class SlurmRealAdapterPreconditionError(SlurmPolicyError):
    """Raised when real-policy mode preconditions are not satisfied."""


def _is_non_empty_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _load_policy(path: Path) -> Mapping[str, Any]:
    if not path.exists():
        raise SlurmPolicyConfigError(f"slurm policy file does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise SlurmPolicyConfigError(
            f"failed to read slurm policy file: {path}"
        ) from exc
    except json.JSONDecodeError as exc:
        raise SlurmPolicyConfigError(f"invalid JSON in slurm policy file: {path}") from exc
    if not isinstance(payload, dict):
        raise SlurmPolicyConfigError("slurm policy top-level object must be a JSON object")
    return cast(Mapping[str, Any], payload)


def parse_slurm_adapter_mode(value: str | None) -> SlurmAdapterMode:
    normalized = (value or "stub-policy").strip().lower()
    if normalized in {"passthrough", "stub-policy", "real-policy"}:
        return cast(SlurmAdapterMode, normalized)
    raise SlurmPolicyConfigError(
        "slurm adapter mode must be one of: passthrough, stub-policy, real-policy"
    )


def parse_slurm_adapter_rollback_guard(value: str | None) -> SlurmAdapterRollbackGuard:
    normalized = (value or "allow").strip().lower()
    if normalized in {"allow", "force-stub-policy", "force-passthrough"}:
        return cast(SlurmAdapterRollbackGuard, normalized)
    raise SlurmPolicyConfigError(
        "slurm adapter rollback guard must be one of: allow, force-stub-policy, force-passthrough"
    )


def resolve_effective_slurm_adapter_mode(
    configured_mode: SlurmAdapterMode,
    *,
    rollback_guard: SlurmAdapterRollbackGuard,
) -> SlurmAdapterMode:
    if rollback_guard == "force-stub-policy":
        return "stub-policy"
    if rollback_guard == "force-passthrough":
        return "passthrough"
    return configured_mode


def _parse_queue_mappings(policy: Mapping[str, Any]) -> dict[str, SlurmQueueMapping]:
    raw_mappings = policy.get("queue_mappings")
    if not isinstance(raw_mappings, list) or not raw_mappings:
        raise SlurmPolicyConfigError("slurm policy queue_mappings must be a non-empty list")

    mappings: dict[str, SlurmQueueMapping] = {}
    for index, raw_mapping in enumerate(raw_mappings):
        if not isinstance(raw_mapping, dict):
            raise SlurmPolicyConfigError(f"queue_mappings[{index}] must be an object")
        queue = raw_mapping.get("queue")
        partition = raw_mapping.get("partition")
        account = raw_mapping.get("account")
        qos = raw_mapping.get("qos")
        max_walltime_minutes = raw_mapping.get("max_walltime_minutes")

        if not _is_non_empty_text(queue):
            raise SlurmPolicyConfigError(
                f"queue_mappings[{index}].queue must be a non-empty string"
            )
        if not _is_non_empty_text(partition):
            raise SlurmPolicyConfigError(
                f"queue_mappings[{index}].partition must be a non-empty string"
            )
        if not _is_non_empty_text(account):
            raise SlurmPolicyConfigError(
                f"queue_mappings[{index}].account must be a non-empty string"
            )
        if not _is_non_empty_text(qos):
            raise SlurmPolicyConfigError(
                f"queue_mappings[{index}].qos must be a non-empty string"
            )
        if not isinstance(max_walltime_minutes, int) or max_walltime_minutes <= 0:
            raise SlurmPolicyConfigError(
                f"queue_mappings[{index}].max_walltime_minutes must be a positive integer"
            )

        queue_name = cast(str, queue).strip()
        if queue_name in mappings:
            raise SlurmPolicyConfigError(f"duplicate queue mapping found for '{queue_name}'")
        mappings[queue_name] = SlurmQueueMapping(
            queue=queue_name,
            partition=cast(str, partition).strip(),
            account=cast(str, account).strip(),
            qos=cast(str, qos).strip(),
            max_walltime_minutes=max_walltime_minutes,
        )
    return mappings


def _resolve_unknown_queue(
    requested_queue: str,
    *,
    policy: Mapping[str, Any],
    mappings: Mapping[str, SlurmQueueMapping],
) -> SlurmQueueResolution:
    mode, default_queue_name = _validate_fallback_policy(policy=policy, mappings=mappings)
    if mode == "deny":
        raise SlurmPolicyDeniedError(
            f"requested queue '{requested_queue}' is not allowed by slurm fallback policy"
        )

    mapping = mappings[default_queue_name]
    return SlurmQueueResolution(
        requested_queue=requested_queue,
        resolved_queue=default_queue_name,
        used_fallback=True,
        mapping=mapping,
    )


def _validate_fallback_policy(
    *, policy: Mapping[str, Any], mappings: Mapping[str, SlurmQueueMapping]
) -> tuple[Literal["deny", "route-default"], str]:
    fallback_policy = policy.get("fallback_policy")
    if not isinstance(fallback_policy, dict):
        raise SlurmPolicyConfigError("fallback_policy must be an object")

    mode = fallback_policy.get("mode")
    if mode == "route-default":
        default_queue = fallback_policy.get("default_queue")
        if not _is_non_empty_text(default_queue):
            raise SlurmPolicyConfigError(
                "fallback_policy.default_queue must be set for route-default mode"
            )
        default_queue_name = cast(str, default_queue).strip()
        mapping = mappings.get(default_queue_name)
        if mapping is None:
            raise SlurmPolicyConfigError(
                "fallback_policy.default_queue must match a queue_mappings.queue"
            )
        return ("route-default", default_queue_name)

    if mode == "deny":
        return ("deny", "")

    raise SlurmPolicyConfigError(
        "fallback_policy.mode must be one of: deny, route-default"
    )


def resolve_runtime_slurm_queue(
    requested_queue: str,
    *,
    policy_path: Path,
) -> SlurmQueueResolution:
    queue_name = requested_queue.strip()
    if not queue_name:
        raise SlurmPolicyDeniedError("requested queue must be a non-empty string")

    policy = _load_policy(policy_path)
    mappings = _parse_queue_mappings(policy)
    mapped = mappings.get(queue_name)
    if mapped is not None:
        return SlurmQueueResolution(
            requested_queue=queue_name,
            resolved_queue=queue_name,
            used_fallback=False,
            mapping=mapped,
        )

    return _resolve_unknown_queue(queue_name, policy=policy, mappings=mappings)


def validate_slurm_policy_file(*, policy_path: Path) -> None:
    policy = _load_policy(policy_path)
    mappings = _parse_queue_mappings(policy)
    _validate_fallback_policy(policy=policy, mappings=mappings)


def resolve_slurm_adapter_stub(
    requested_queue: str,
    *,
    policy_path: Path | None,
) -> SlurmAdapterBoundaryStub:
    queue_name = requested_queue.strip()
    if not queue_name:
        raise SlurmPolicyDeniedError("requested queue must be a non-empty string")

    if policy_path is None:
        return SlurmAdapterBoundaryStub(
            adapter="passthrough",
            contract_version=SLURM_ADAPTER_CONTRACT_VERSION,
            requested_queue=queue_name,
            resolved_queue=queue_name,
            used_fallback=False,
            mapping=None,
        )

    resolution = resolve_runtime_slurm_queue(queue_name, policy_path=policy_path)
    return SlurmAdapterBoundaryStub(
        adapter="stub-policy",
        contract_version=SLURM_ADAPTER_CONTRACT_VERSION,
        requested_queue=resolution.requested_queue,
        resolved_queue=resolution.resolved_queue,
        used_fallback=resolution.used_fallback,
        mapping=resolution.mapping,
    )


def resolve_slurm_adapter_real(
    requested_queue: str,
    *,
    policy_path: Path | None,
    command_timeout_seconds: int = 5,
) -> SlurmAdapterBoundaryStub:
    queue_name = requested_queue.strip()
    if not queue_name:
        raise SlurmPolicyDeniedError("requested queue must be a non-empty string")
    if policy_path is None:
        raise SlurmRealAdapterPreconditionError(
            "real-policy requires ZPE_SLURM_POLICY_PATH to be configured"
        )

    validate_slurm_policy_file(policy_path=policy_path)
    try:
        probe = subprocess.run(
            ["scontrol", "ping"],
            capture_output=True,
            text=True,
            check=False,
            timeout=command_timeout_seconds,
        )
    except FileNotFoundError as exc:
        raise SlurmRealAdapterPreconditionError(
            "real-policy requires scontrol command to be available"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise SlurmRealAdapterPreconditionError(
            f"real-policy precondition probe timed out after {command_timeout_seconds}s: scontrol ping"
        ) from exc
    except OSError as exc:
        raise SlurmRealAdapterPreconditionError(
            f"real-policy precondition probe failed to launch: {exc}"
        ) from exc

    if probe.returncode != 0:
        detail = (probe.stderr or probe.stdout or "").strip()
        if detail:
            raise SlurmRealAdapterPreconditionError(
                f"real-policy precondition probe failed: {detail}"
            )
        raise SlurmRealAdapterPreconditionError(
            "real-policy precondition probe failed: scontrol ping returned non-zero"
        )

    resolution = resolve_runtime_slurm_queue(queue_name, policy_path=policy_path)
    return SlurmAdapterBoundaryStub(
        adapter="real-policy",
        contract_version=SLURM_ADAPTER_CONTRACT_VERSION,
        requested_queue=resolution.requested_queue,
        resolved_queue=resolution.resolved_queue,
        used_fallback=resolution.used_fallback,
        mapping=resolution.mapping,
    )
