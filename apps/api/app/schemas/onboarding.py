from __future__ import annotations

from typing import Any, Literal

from pydantic import Field

from .base import ApiModel


OnboardingSchemaVersion = Literal["byo_onboarding_manifest/v1"]
OnboardingPathMode = Literal[
    "path-a-existing-slurm",
    "path-b-bootstrap-slurm-aiida",
]
OnboardingIssueSection = Literal[
    "manifest",
    "policy",
    "registration",
    "queue_resolution",
]


class OnboardingDryRunRequest(ApiModel):
    schema_version: OnboardingSchemaVersion = "byo_onboarding_manifest/v1"
    path_mode: OnboardingPathMode
    policy: dict[str, Any]
    registration: dict[str, Any] | None = None
    requested_queue: str | None = None


class OnboardingValidationIssue(ApiModel):
    section: OnboardingIssueSection
    field_path: str | None = None
    code: str
    message: str
    action: str


class OnboardingQueueDecision(ApiModel):
    requested_queue: str
    resolved_queue: str
    used_fallback: bool


class OnboardingDryRunReport(ApiModel):
    policy_valid: bool
    registration_valid: bool
    queue_resolution_valid: bool
    error_count: int
    policy_error_count: int
    registration_error_count: int
    queue_resolution_error_count: int
    queue_decision: OnboardingQueueDecision | None = None


class OnboardingDryRunResponse(ApiModel):
    status: Literal["valid", "invalid"]
    normalized_manifest: OnboardingDryRunRequest
    report: OnboardingDryRunReport
    errors: list[OnboardingValidationIssue] = Field(default_factory=list)
