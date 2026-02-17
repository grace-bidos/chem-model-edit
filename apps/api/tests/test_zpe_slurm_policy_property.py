from __future__ import annotations

import json
from pathlib import Path
import tempfile

from hypothesis import assume, given
from hypothesis import strategies as st
import pytest

from services.zpe.slurm_policy import SlurmPolicyDeniedError, resolve_runtime_slurm_queue

_QUEUE_RE = r"[a-z][a-z0-9_-]{0,10}"


def _write_policy(
    path: Path,
    *,
    queue_names: list[str],
    fallback_mode: str,
    default_queue: str,
) -> None:
    queue_mappings = [
        {
            "queue": queue,
            "partition": f"partition-{idx}",
            "account": f"account-{idx}",
            "qos": f"qos-{idx}",
            "max_walltime_minutes": 30 + idx,
        }
        for idx, queue in enumerate(queue_names)
    ]
    fallback_policy: dict[str, str] = {"mode": fallback_mode}
    if fallback_mode == "route-default":
        fallback_policy["default_queue"] = default_queue
    path.write_text(
        json.dumps(
            {
                "queue_mappings": queue_mappings,
                "fallback_policy": fallback_policy,
            }
        ),
        encoding="utf-8",
    )


@given(
    queue_names=st.lists(
        st.from_regex(_QUEUE_RE, fullmatch=True),
        min_size=1,
        max_size=8,
        unique=True,
    ),
    requested_index=st.integers(min_value=0, max_value=7),
)
def test_resolve_runtime_slurm_queue_keeps_known_queue_mapping(
    queue_names: list[str], requested_index: int
) -> None:
    assume(requested_index < len(queue_names))
    requested_queue = queue_names[requested_index]
    with tempfile.TemporaryDirectory() as tmp_dir:
        policy_path = Path(tmp_dir) / "policy.json"
        _write_policy(
            policy_path,
            queue_names=queue_names,
            fallback_mode="route-default",
            default_queue=queue_names[0],
        )
        resolution = resolve_runtime_slurm_queue(requested_queue, policy_path=policy_path)

    assert resolution.requested_queue == requested_queue
    assert resolution.resolved_queue == requested_queue
    assert resolution.used_fallback is False
    assert resolution.mapping.queue == requested_queue


@given(
    queue_names=st.lists(
        st.from_regex(_QUEUE_RE, fullmatch=True),
        min_size=1,
        max_size=8,
        unique=True,
    ),
    requested_queue=st.from_regex(_QUEUE_RE, fullmatch=True),
    fallback_mode=st.sampled_from(["route-default", "deny"]),
)
def test_resolve_runtime_slurm_queue_unknown_queue_obeys_fallback_mode(
    queue_names: list[str],
    requested_queue: str,
    fallback_mode: str,
) -> None:
    assume(requested_queue not in queue_names)
    default_queue = queue_names[0]
    with tempfile.TemporaryDirectory() as tmp_dir:
        policy_path = Path(tmp_dir) / "policy.json"
        _write_policy(
            policy_path,
            queue_names=queue_names,
            fallback_mode=fallback_mode,
            default_queue=default_queue,
        )

        if fallback_mode == "deny":
            with pytest.raises(SlurmPolicyDeniedError):
                resolve_runtime_slurm_queue(requested_queue, policy_path=policy_path)
            return

        resolution = resolve_runtime_slurm_queue(requested_queue, policy_path=policy_path)
    assert resolution.requested_queue == requested_queue
    assert resolution.resolved_queue == default_queue
    assert resolution.used_fallback is True
    assert resolution.mapping.queue == default_queue
