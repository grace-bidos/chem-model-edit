from __future__ import annotations

import json
from pathlib import Path

import pytest

from services.zpe.slurm_policy import (
    SLURM_ADAPTER_CONTRACT_VERSION,
    SlurmPolicyConfigError,
    SlurmPolicyDeniedError,
    resolve_slurm_adapter_stub,
    resolve_runtime_slurm_queue,
)


def _write_policy(path: Path, *, fallback_mode: str = "route-default") -> None:
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
            },
            {
                "queue": "highmem",
                "partition": "long",
                "account": "chem-premium",
                "qos": "priority",
                "max_walltime_minutes": 720,
            },
        ],
        "fallback_policy": {"mode": fallback_mode, "default_queue": "standard"},
    }
    if fallback_mode == "deny":
        policy["fallback_policy"] = {"mode": "deny"}
    path.write_text(json.dumps(policy), encoding="utf-8")


def test_resolve_runtime_slurm_queue_known_queue(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    _write_policy(policy_path)

    resolution = resolve_runtime_slurm_queue("highmem", policy_path=policy_path)

    assert resolution.requested_queue == "highmem"
    assert resolution.resolved_queue == "highmem"
    assert resolution.used_fallback is False
    assert resolution.mapping.partition == "long"
    assert resolution.mapping.account == "chem-premium"
    assert resolution.mapping.qos == "priority"


def test_resolve_runtime_slurm_queue_unknown_queue_route_default(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    _write_policy(policy_path, fallback_mode="route-default")

    resolution = resolve_runtime_slurm_queue("legacy", policy_path=policy_path)

    assert resolution.requested_queue == "legacy"
    assert resolution.resolved_queue == "standard"
    assert resolution.used_fallback is True
    assert resolution.mapping.partition == "short"
    assert resolution.mapping.account == "chem-default"
    assert resolution.mapping.qos == "normal"


def test_resolve_runtime_slurm_queue_unknown_queue_deny(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    _write_policy(policy_path, fallback_mode="deny")

    with pytest.raises(SlurmPolicyDeniedError):
        resolve_runtime_slurm_queue("legacy", policy_path=policy_path)


def test_resolve_runtime_slurm_queue_wraps_policy_read_errors(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy-dir"
    policy_path.mkdir()

    with pytest.raises(SlurmPolicyConfigError):
        resolve_runtime_slurm_queue("standard", policy_path=policy_path)


def test_resolve_slurm_adapter_stub_passthrough_without_policy() -> None:
    resolution = resolve_slurm_adapter_stub(" standard ", policy_path=None)

    assert resolution.adapter == "passthrough"
    assert resolution.contract_version == SLURM_ADAPTER_CONTRACT_VERSION
    assert resolution.requested_queue == "standard"
    assert resolution.resolved_queue == "standard"
    assert resolution.used_fallback is False
    assert resolution.mapping is None


def test_resolve_slurm_adapter_stub_uses_policy_when_present(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    _write_policy(policy_path, fallback_mode="route-default")

    resolution = resolve_slurm_adapter_stub("legacy", policy_path=policy_path)

    assert resolution.adapter == "stub-policy"
    assert resolution.contract_version == SLURM_ADAPTER_CONTRACT_VERSION
    assert resolution.requested_queue == "legacy"
    assert resolution.resolved_queue == "standard"
    assert resolution.used_fallback is True
    assert resolution.mapping is not None
    assert resolution.mapping.partition == "short"
