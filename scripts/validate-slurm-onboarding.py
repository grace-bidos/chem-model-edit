#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ALLOWED_FALLBACK_MODES = {"deny", "route-default"}
BLOCKED_STATES = {"down", "drain", "draining", "fail"}


def _is_non_empty_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ValueError(f"file does not exist: {path}")
    try:
        content = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in {path}: {exc}") from exc
    if not isinstance(content, dict):
        raise ValueError(f"expected object at top-level in {path}")
    return content


def _validate_text_list(data: Any, field_name: str, errors: list[str]) -> list[str]:
    if not isinstance(data, list) or not data:
        errors.append(f"{field_name} must be a non-empty list")
        return []
    values: list[str] = []
    for index, item in enumerate(data):
        if not _is_non_empty_text(item):
            errors.append(f"{field_name}[{index}] must be a non-empty string")
            continue
        values.append(item.strip())
    if len(set(values)) != len(values):
        errors.append(f"{field_name} must not contain duplicates")
    return values


def validate_policy(policy: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    if not _is_non_empty_text(policy.get("cluster_name")):
        errors.append("cluster_name must be a non-empty string")

    slurm = policy.get("slurm")
    if not isinstance(slurm, dict):
        errors.append("slurm must be an object")
        slurm = {}

    partitions = _validate_text_list(slurm.get("partitions"), "slurm.partitions", errors)
    accounts = _validate_text_list(slurm.get("accounts"), "slurm.accounts", errors)
    qos_values = _validate_text_list(slurm.get("qos"), "slurm.qos", errors)

    queue_mappings = policy.get("queue_mappings")
    if not isinstance(queue_mappings, list) or not queue_mappings:
        errors.append("queue_mappings must be a non-empty list")
        queue_mappings = []

    queue_names: list[str] = []
    for index, mapping in enumerate(queue_mappings):
        prefix = f"queue_mappings[{index}]"
        if not isinstance(mapping, dict):
            errors.append(f"{prefix} must be an object")
            continue

        queue = mapping.get("queue")
        partition = mapping.get("partition")
        account = mapping.get("account")
        qos = mapping.get("qos")
        max_walltime_minutes = mapping.get("max_walltime_minutes")

        if not _is_non_empty_text(queue):
            errors.append(f"{prefix}.queue must be a non-empty string")
        else:
            queue_names.append(queue.strip())

        if not _is_non_empty_text(partition):
            errors.append(f"{prefix}.partition must be a non-empty string")
        elif partitions and partition not in partitions:
            errors.append(f"{prefix}.partition must be one of slurm.partitions")

        if not _is_non_empty_text(account):
            errors.append(f"{prefix}.account must be a non-empty string")
        elif accounts and account not in accounts:
            errors.append(f"{prefix}.account must be one of slurm.accounts")

        if not _is_non_empty_text(qos):
            errors.append(f"{prefix}.qos must be a non-empty string")
        elif qos_values and qos not in qos_values:
            errors.append(f"{prefix}.qos must be one of slurm.qos")

        if not isinstance(max_walltime_minutes, int) or max_walltime_minutes <= 0:
            errors.append(f"{prefix}.max_walltime_minutes must be a positive integer")

    if len(set(queue_names)) != len(queue_names):
        errors.append("queue_mappings queue names must be unique")

    fallback_policy = policy.get("fallback_policy")
    if not isinstance(fallback_policy, dict):
        errors.append("fallback_policy must be an object")
        fallback_policy = {}

    mode = fallback_policy.get("mode")
    if mode not in ALLOWED_FALLBACK_MODES:
        errors.append(
            "fallback_policy.mode must be one of: "
            + ", ".join(sorted(ALLOWED_FALLBACK_MODES))
        )
    if mode == "route-default":
        default_queue = fallback_policy.get("default_queue")
        if not _is_non_empty_text(default_queue):
            errors.append("fallback_policy.default_queue must be set for route-default mode")
        elif default_queue not in queue_names:
            errors.append("fallback_policy.default_queue must match a queue_mappings.queue")

    registration_requirements = policy.get("registration_requirements")
    if not isinstance(registration_requirements, dict):
        errors.append("registration_requirements must be an object")
        registration_requirements = {}

    _validate_text_list(
        registration_requirements.get("required_node_labels"),
        "registration_requirements.required_node_labels",
        errors,
    )
    _validate_text_list(
        registration_requirements.get("required_health_checks"),
        "registration_requirements.required_health_checks",
        errors,
    )

    deny_if_drain_or_down = registration_requirements.get("deny_if_drain_or_down")
    if not isinstance(deny_if_drain_or_down, bool):
        errors.append("registration_requirements.deny_if_drain_or_down must be boolean")

    return errors


def validate_registration(
    policy: dict[str, Any], registration: dict[str, Any]
) -> list[str]:
    errors: list[str] = []

    if registration.get("cluster_name") != policy.get("cluster_name"):
        errors.append("registration.cluster_name must match policy.cluster_name")

    enabled_queues = registration.get("enabled_queues")
    if not isinstance(enabled_queues, list) or not enabled_queues:
        errors.append("registration.enabled_queues must be a non-empty list")
        enabled_queues = []
    else:
        for index, queue in enumerate(enabled_queues):
            if not _is_non_empty_text(queue):
                errors.append(f"registration.enabled_queues[{index}] must be a non-empty string")

    mapped_queues = {
        mapping.get("queue")
        for mapping in policy.get("queue_mappings", [])
        if isinstance(mapping, dict) and _is_non_empty_text(mapping.get("queue"))
    }
    for queue in enabled_queues:
        if isinstance(queue, str) and queue not in mapped_queues:
            errors.append(f"registration.enabled_queues includes unmapped queue: {queue}")

    requirements = policy.get("registration_requirements", {})
    required_labels = requirements.get("required_node_labels", [])
    labels = registration.get("node_labels")
    if not isinstance(labels, dict):
        errors.append("registration.node_labels must be an object")
        labels = {}
    for label in required_labels:
        value = labels.get(label)
        if not _is_non_empty_text(value):
            errors.append(f"registration.node_labels.{label} must be a non-empty string")

    required_checks = requirements.get("required_health_checks", [])
    checks = registration.get("health_checks")
    if not isinstance(checks, dict):
        errors.append("registration.health_checks must be an object")
        checks = {}
    for check_name in required_checks:
        if checks.get(check_name) is not True:
            errors.append(f"registration.health_checks.{check_name} must be true")

    deny_if_drain_or_down = requirements.get("deny_if_drain_or_down", False)
    slurm_state = registration.get("slurm_state")
    if not _is_non_empty_text(slurm_state):
        errors.append("registration.slurm_state must be a non-empty string")
    elif deny_if_drain_or_down and slurm_state.strip().lower() in BLOCKED_STATES:
        errors.append(
            "registration.slurm_state is not allowed when "
            "deny_if_drain_or_down=true"
        )

    return errors


def resolve_queue(policy: dict[str, Any], requested_queue: str) -> tuple[str, bool]:
    queue_names = {
        mapping.get("queue")
        for mapping in policy.get("queue_mappings", [])
        if isinstance(mapping, dict) and _is_non_empty_text(mapping.get("queue"))
    }
    if requested_queue in queue_names:
        return requested_queue, False

    fallback_policy = policy.get("fallback_policy", {})
    mode = fallback_policy.get("mode")
    if mode == "route-default":
        default_queue = fallback_policy.get("default_queue")
        if _is_non_empty_text(default_queue):
            return default_queue, True
    raise ValueError(
        f"requested queue '{requested_queue}' is unknown and fallback mode does not allow routing"
    )


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Slurm BYO onboarding policy and optional node registration payload."
        )
    )
    parser.add_argument(
        "--policy",
        required=True,
        help="Path to onboarding policy JSON file",
    )
    parser.add_argument(
        "--registration",
        help="Path to node registration JSON file",
    )
    parser.add_argument(
        "--requested-queue",
        help="Optional queue name to test fallback behavior",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    all_errors: list[str] = []

    try:
        policy = _load_json(Path(args.policy))
    except ValueError as exc:
        print(f"ERROR: {exc}")
        return 1

    policy_errors = validate_policy(policy)
    all_errors.extend([f"policy: {err}" for err in policy_errors])

    registration: dict[str, Any] | None = None
    if args.registration:
        try:
            registration = _load_json(Path(args.registration))
        except ValueError as exc:
            all_errors.append(f"registration: {exc}")

    if registration is not None:
        registration_errors = validate_registration(policy, registration)
        all_errors.extend([f"registration: {err}" for err in registration_errors])

    if args.requested_queue:
        try:
            resolved_queue, used_fallback = resolve_queue(policy, args.requested_queue)
            print(
                "QUEUE_DECISION "
                f"requested={args.requested_queue} "
                f"resolved={resolved_queue} "
                f"used_fallback={'yes' if used_fallback else 'no'}"
            )
        except ValueError as exc:
            all_errors.append(f"queue_resolution: {exc}")

    if all_errors:
        for error in all_errors:
            print(f"ERROR: {error}")
        print(f"Validation failed with {len(all_errors)} error(s).")
        return 1

    print("Validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
