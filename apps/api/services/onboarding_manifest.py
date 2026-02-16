from __future__ import annotations

from typing import Any

from app.schemas.onboarding import (
    OnboardingDryRunReport,
    OnboardingDryRunRequest,
    OnboardingDryRunResponse,
    OnboardingIssueSection,
    OnboardingQueueDecision,
    OnboardingValidationIssue,
)

ALLOWED_FALLBACK_MODES = {"deny", "route-default"}
BLOCKED_STATES = {"down", "drain", "draining", "fail"}


def _is_non_empty_text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _normalize_json(value: Any) -> Any:
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, list):
        return [_normalize_json(item) for item in value]
    if isinstance(value, dict):
        normalized: dict[str, Any] = {}
        for key, item in value.items():
            normalized[key] = _normalize_json(item)
        return normalized
    return value


def _issue(
    *,
    section: OnboardingIssueSection,
    field_path: str | None,
    message: str,
    action: str,
    code: str = "invalid_value",
) -> OnboardingValidationIssue:
    return OnboardingValidationIssue(
        section=section,
        field_path=field_path,
        code=code,
        message=message,
        action=action,
    )


def _validate_text_list(
    *,
    section: OnboardingIssueSection,
    data: Any,
    field_name: str,
    errors: list[OnboardingValidationIssue],
) -> list[str]:
    if not isinstance(data, list) or not data:
        errors.append(
            _issue(
                section=section,
                field_path=field_name,
                message=f"{field_name} must be a non-empty list",
                action="Provide at least one non-empty string value.",
                code="missing_or_invalid_list",
            )
        )
        return []

    values: list[str] = []
    for index, item in enumerate(data):
        if not _is_non_empty_text(item):
            errors.append(
                _issue(
                    section=section,
                    field_path=f"{field_name}[{index}]",
                    message=f"{field_name}[{index}] must be a non-empty string",
                    action="Replace this entry with a non-empty string.",
                    code="invalid_item_type",
                )
            )
            continue
        values.append(item.strip())

    if len(set(values)) != len(values):
        errors.append(
            _issue(
                section=section,
                field_path=field_name,
                message=f"{field_name} must not contain duplicates",
                action="Remove duplicated values from this list.",
                code="duplicate_items",
            )
        )
    return values


def _validate_policy(policy: dict[str, Any]) -> list[OnboardingValidationIssue]:
    errors: list[OnboardingValidationIssue] = []

    if not _is_non_empty_text(policy.get("cluster_name")):
        errors.append(
            _issue(
                section="policy",
                field_path="cluster_name",
                message="cluster_name must be a non-empty string",
                action="Set cluster_name to the Slurm cluster identifier.",
            )
        )

    slurm = policy.get("slurm")
    if not isinstance(slurm, dict):
        errors.append(
            _issue(
                section="policy",
                field_path="slurm",
                message="slurm must be an object",
                action="Provide slurm with partitions/accounts/qos lists.",
                code="invalid_object",
            )
        )
        slurm = {}

    partitions = _validate_text_list(
        section="policy",
        data=slurm.get("partitions"),
        field_name="slurm.partitions",
        errors=errors,
    )
    accounts = _validate_text_list(
        section="policy",
        data=slurm.get("accounts"),
        field_name="slurm.accounts",
        errors=errors,
    )
    qos_values = _validate_text_list(
        section="policy",
        data=slurm.get("qos"),
        field_name="slurm.qos",
        errors=errors,
    )

    queue_mappings = policy.get("queue_mappings")
    if not isinstance(queue_mappings, list) or not queue_mappings:
        errors.append(
            _issue(
                section="policy",
                field_path="queue_mappings",
                message="queue_mappings must be a non-empty list",
                action="Add at least one queue mapping entry.",
                code="missing_or_invalid_list",
            )
        )
        queue_mappings = []

    queue_names: list[str] = []
    for index, mapping in enumerate(queue_mappings):
        prefix = f"queue_mappings[{index}]"
        if not isinstance(mapping, dict):
            errors.append(
                _issue(
                    section="policy",
                    field_path=prefix,
                    message=f"{prefix} must be an object",
                    action="Provide queue, partition, account, qos, and walltime fields.",
                    code="invalid_object",
                )
            )
            continue

        queue = mapping.get("queue")
        partition = mapping.get("partition")
        account = mapping.get("account")
        qos = mapping.get("qos")
        max_walltime_minutes = mapping.get("max_walltime_minutes")

        if not _is_non_empty_text(queue):
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.queue",
                    message=f"{prefix}.queue must be a non-empty string",
                    action="Set a stable logical queue name.",
                )
            )
        else:
            assert isinstance(queue, str)
            queue_names.append(queue.strip())

        if not _is_non_empty_text(partition):
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.partition",
                    message=f"{prefix}.partition must be a non-empty string",
                    action="Set a Slurm partition name.",
                )
            )
        elif partitions and partition not in partitions:
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.partition",
                    message=f"{prefix}.partition must be one of slurm.partitions",
                    action="Use a partition from slurm.partitions.",
                    code="allowlist_mismatch",
                )
            )

        if not _is_non_empty_text(account):
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.account",
                    message=f"{prefix}.account must be a non-empty string",
                    action="Set a Slurm account name.",
                )
            )
        elif accounts and account not in accounts:
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.account",
                    message=f"{prefix}.account must be one of slurm.accounts",
                    action="Use an account from slurm.accounts.",
                    code="allowlist_mismatch",
                )
            )

        if not _is_non_empty_text(qos):
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.qos",
                    message=f"{prefix}.qos must be a non-empty string",
                    action="Set a Slurm QoS name.",
                )
            )
        elif qos_values and qos not in qos_values:
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.qos",
                    message=f"{prefix}.qos must be one of slurm.qos",
                    action="Use a QoS from slurm.qos.",
                    code="allowlist_mismatch",
                )
            )

        if (
            isinstance(max_walltime_minutes, bool)
            or not isinstance(max_walltime_minutes, int)
            or max_walltime_minutes <= 0
        ):
            errors.append(
                _issue(
                    section="policy",
                    field_path=f"{prefix}.max_walltime_minutes",
                    message=f"{prefix}.max_walltime_minutes must be a positive integer",
                    action="Set max_walltime_minutes to an integer greater than 0.",
                    code="invalid_integer",
                )
            )

    if len(set(queue_names)) != len(queue_names):
        errors.append(
            _issue(
                section="policy",
                field_path="queue_mappings",
                message="queue_mappings queue names must be unique",
                action="Ensure each queue_mappings[].queue appears only once.",
                code="duplicate_items",
            )
        )

    fallback_policy = policy.get("fallback_policy")
    if not isinstance(fallback_policy, dict):
        errors.append(
            _issue(
                section="policy",
                field_path="fallback_policy",
                message="fallback_policy must be an object",
                action="Provide mode and optional default_queue values.",
                code="invalid_object",
            )
        )
        fallback_policy = {}

    mode = fallback_policy.get("mode")
    if mode not in ALLOWED_FALLBACK_MODES:
        errors.append(
            _issue(
                section="policy",
                field_path="fallback_policy.mode",
                message=(
                    "fallback_policy.mode must be one of: "
                    + ", ".join(sorted(ALLOWED_FALLBACK_MODES))
                ),
                action="Use either deny or route-default.",
                code="invalid_choice",
            )
        )
    if mode == "route-default":
        default_queue = fallback_policy.get("default_queue")
        if not _is_non_empty_text(default_queue):
            errors.append(
                _issue(
                    section="policy",
                    field_path="fallback_policy.default_queue",
                    message=(
                        "fallback_policy.default_queue must be set for route-default mode"
                    ),
                    action="Set default_queue to an existing mapped queue.",
                )
            )
        elif default_queue not in queue_names:
            errors.append(
                _issue(
                    section="policy",
                    field_path="fallback_policy.default_queue",
                    message="fallback_policy.default_queue must match a queue_mappings.queue",
                    action="Choose default_queue from queue_mappings[].queue.",
                    code="allowlist_mismatch",
                )
            )

    registration_requirements = policy.get("registration_requirements")
    if not isinstance(registration_requirements, dict):
        errors.append(
            _issue(
                section="policy",
                field_path="registration_requirements",
                message="registration_requirements must be an object",
                action="Provide required labels/checks and deny_if_drain_or_down.",
                code="invalid_object",
            )
        )
        registration_requirements = {}

    _validate_text_list(
        section="policy",
        data=registration_requirements.get("required_node_labels"),
        field_name="registration_requirements.required_node_labels",
        errors=errors,
    )
    _validate_text_list(
        section="policy",
        data=registration_requirements.get("required_health_checks"),
        field_name="registration_requirements.required_health_checks",
        errors=errors,
    )

    deny_if_drain_or_down = registration_requirements.get("deny_if_drain_or_down")
    if not isinstance(deny_if_drain_or_down, bool):
        errors.append(
            _issue(
                section="policy",
                field_path="registration_requirements.deny_if_drain_or_down",
                message="registration_requirements.deny_if_drain_or_down must be boolean",
                action="Set deny_if_drain_or_down to true or false.",
                code="invalid_boolean",
            )
        )

    return errors


def _validate_registration(
    policy: dict[str, Any], registration: dict[str, Any]
) -> list[OnboardingValidationIssue]:
    errors: list[OnboardingValidationIssue] = []

    if registration.get("cluster_name") != policy.get("cluster_name"):
        errors.append(
            _issue(
                section="registration",
                field_path="cluster_name",
                message="registration.cluster_name must match policy.cluster_name",
                action="Set registration.cluster_name to the same value as policy.cluster_name.",
                code="mismatch",
            )
        )

    enabled_queues = registration.get("enabled_queues")
    if not isinstance(enabled_queues, list) or not enabled_queues:
        errors.append(
            _issue(
                section="registration",
                field_path="enabled_queues",
                message="registration.enabled_queues must be a non-empty list",
                action="Add one or more mapped queue names.",
                code="missing_or_invalid_list",
            )
        )
        enabled_queues = []
    else:
        for index, queue in enumerate(enabled_queues):
            if not _is_non_empty_text(queue):
                errors.append(
                    _issue(
                        section="registration",
                        field_path=f"enabled_queues[{index}]",
                        message=(
                            "registration.enabled_queues"
                            f"[{index}] must be a non-empty string"
                        ),
                        action="Replace this entry with a non-empty mapped queue name.",
                    )
                )

    mapped_queues = {
        mapping.get("queue")
        for mapping in policy.get("queue_mappings", [])
        if isinstance(mapping, dict) and _is_non_empty_text(mapping.get("queue"))
    }

    for queue in enabled_queues:
        if isinstance(queue, str) and queue not in mapped_queues:
            errors.append(
                _issue(
                    section="registration",
                    field_path="enabled_queues",
                    message=(
                        "registration.enabled_queues includes unmapped queue: "
                        + queue
                    ),
                    action="Use queue names defined in policy.queue_mappings.",
                    code="allowlist_mismatch",
                )
            )

    requirements = policy.get("registration_requirements", {})
    required_labels = requirements.get("required_node_labels", [])
    labels = registration.get("node_labels")
    if not isinstance(labels, dict):
        errors.append(
            _issue(
                section="registration",
                field_path="node_labels",
                message="registration.node_labels must be an object",
                action="Provide key/value labels for required_node_labels.",
                code="invalid_object",
            )
        )
        labels = {}
    for label in required_labels if isinstance(required_labels, list) else []:
        value = labels.get(label)
        if not _is_non_empty_text(value):
            errors.append(
                _issue(
                    section="registration",
                    field_path=f"node_labels.{label}",
                    message=(
                        "registration.node_labels."
                        + str(label)
                        + " must be a non-empty string"
                    ),
                    action="Set this required node label to a non-empty string.",
                )
            )

    required_checks = requirements.get("required_health_checks", [])
    checks = registration.get("health_checks")
    if not isinstance(checks, dict):
        errors.append(
            _issue(
                section="registration",
                field_path="health_checks",
                message="registration.health_checks must be an object",
                action="Provide boolean values for required health checks.",
                code="invalid_object",
            )
        )
        checks = {}
    for check_name in required_checks if isinstance(required_checks, list) else []:
        if checks.get(check_name) is not True:
            errors.append(
                _issue(
                    section="registration",
                    field_path=f"health_checks.{check_name}",
                    message=(
                        "registration.health_checks."
                        + str(check_name)
                        + " must be true"
                    ),
                    action="Run and pass the required health check before registration.",
                    code="check_failed",
                )
            )

    deny_if_drain_or_down = requirements.get("deny_if_drain_or_down", False)
    slurm_state = registration.get("slurm_state")
    if not _is_non_empty_text(slurm_state):
        errors.append(
            _issue(
                section="registration",
                field_path="slurm_state",
                message="registration.slurm_state must be a non-empty string",
                action="Set current Slurm state (for example, idle).",
            )
        )
    elif (
        deny_if_drain_or_down
        and isinstance(slurm_state, str)
        and slurm_state.strip().lower() in BLOCKED_STATES
    ):
        errors.append(
            _issue(
                section="registration",
                field_path="slurm_state",
                message=(
                    "registration.slurm_state is not allowed when "
                    "deny_if_drain_or_down=true"
                ),
                action="Bring the node to an allowed Slurm state before onboarding.",
                code="blocked_state",
            )
        )

    return errors


def _resolve_queue(policy: dict[str, Any], requested_queue: str) -> tuple[str, bool]:
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
        "requested queue '"
        + requested_queue
        + "' is unknown and fallback mode does not allow routing"
    )


def run_onboarding_dry_run(manifest: OnboardingDryRunRequest) -> OnboardingDryRunResponse:
    normalized_policy = _normalize_json(manifest.policy)
    normalized_registration = _normalize_json(manifest.registration)
    normalized_requested_queue = _normalize_json(manifest.requested_queue)

    normalized_manifest = OnboardingDryRunRequest(
        schema_version=manifest.schema_version,
        path_mode=manifest.path_mode,
        policy=normalized_policy,
        registration=normalized_registration,
        requested_queue=normalized_requested_queue,
    )

    policy_errors = _validate_policy(normalized_policy)
    registration_errors: list[OnboardingValidationIssue] = []
    queue_resolution_errors: list[OnboardingValidationIssue] = []
    queue_decision: OnboardingQueueDecision | None = None

    if normalized_registration is not None:
        registration_errors = _validate_registration(
            normalized_policy,
            normalized_registration,
        )

    if isinstance(normalized_requested_queue, str) and normalized_requested_queue:
        try:
            resolved_queue, used_fallback = _resolve_queue(
                normalized_policy,
                normalized_requested_queue,
            )
            queue_decision = OnboardingQueueDecision(
                requested_queue=normalized_requested_queue,
                resolved_queue=resolved_queue,
                used_fallback=used_fallback,
            )
        except ValueError as exc:
            queue_resolution_errors.append(
                _issue(
                    section="queue_resolution",
                    field_path="requested_queue",
                    message=str(exc),
                    action="Use a mapped queue or switch fallback_policy.mode to route-default.",
                    code="queue_resolution_failed",
                )
            )

    errors = [*policy_errors, *registration_errors, *queue_resolution_errors]
    report = OnboardingDryRunReport(
        policy_valid=not policy_errors,
        registration_valid=not registration_errors,
        queue_resolution_valid=not queue_resolution_errors,
        error_count=len(errors),
        policy_error_count=len(policy_errors),
        registration_error_count=len(registration_errors),
        queue_resolution_error_count=len(queue_resolution_errors),
        queue_decision=queue_decision,
    )
    return OnboardingDryRunResponse(
        status="valid" if not errors else "invalid",
        normalized_manifest=normalized_manifest,
        report=report,
        errors=errors,
    )
