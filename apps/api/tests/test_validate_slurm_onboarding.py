from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = REPO_ROOT / "scripts" / "validate-slurm-onboarding.py"


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def _run_validator(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )


def _base_policy() -> dict:
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
        "fallback_policy": {
            "mode": "route-default",
            "default_queue": "standard",
        },
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


def test_validate_slurm_onboarding_success_with_fallback_resolution(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    registration_path = tmp_path / "registration.json"

    _write_json(policy_path, _base_policy())
    _write_json(registration_path, _base_registration())

    result = _run_validator(
        "--policy",
        str(policy_path),
        "--registration",
        str(registration_path),
        "--requested-queue",
        "experimental",
    )

    assert result.returncode == 0
    assert "Validation passed." in result.stdout
    assert "resolved=standard" in result.stdout
    assert "used_fallback=yes" in result.stdout


def test_validate_slurm_onboarding_deny_policy_blocks_unknown_queue(tmp_path: Path) -> None:
    policy = _base_policy()
    policy["fallback_policy"] = {"mode": "deny"}

    policy_path = tmp_path / "policy.json"
    _write_json(policy_path, policy)

    result = _run_validator(
        "--policy",
        str(policy_path),
        "--requested-queue",
        "experimental",
    )

    assert result.returncode == 1
    assert "queue_resolution:" in result.stdout


def test_validate_slurm_onboarding_rejects_partition_not_in_allowlist(tmp_path: Path) -> None:
    policy = _base_policy()
    policy["queue_mappings"][0]["partition"] = "gpu"

    policy_path = tmp_path / "policy.json"
    _write_json(policy_path, policy)

    result = _run_validator("--policy", str(policy_path))

    assert result.returncode == 1
    assert "queue_mappings[0].partition must be one of slurm.partitions" in result.stdout


def test_validate_slurm_onboarding_rejects_failed_required_health_check(tmp_path: Path) -> None:
    policy_path = tmp_path / "policy.json"
    registration_path = tmp_path / "registration.json"

    registration = _base_registration()
    registration["health_checks"]["sbatch"] = False

    _write_json(policy_path, _base_policy())
    _write_json(registration_path, registration)

    result = _run_validator(
        "--policy",
        str(policy_path),
        "--registration",
        str(registration_path),
    )

    assert result.returncode == 1
    assert "registration.health_checks.sbatch must be true" in result.stdout
