from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, cast


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


class SlurmPolicyError(RuntimeError):
    """Base exception for runtime Slurm policy resolution errors."""


class SlurmPolicyDeniedError(SlurmPolicyError):
    """Raised when policy denies the requested queue."""


class SlurmPolicyConfigError(SlurmPolicyError):
    """Raised when policy content is missing or invalid for runtime resolution."""


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
        return SlurmQueueResolution(
            requested_queue=requested_queue,
            resolved_queue=default_queue_name,
            used_fallback=True,
            mapping=mapping,
        )

    if mode == "deny":
        raise SlurmPolicyDeniedError(
            f"requested queue '{requested_queue}' is not allowed by slurm fallback policy"
        )

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
