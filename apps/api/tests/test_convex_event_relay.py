from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, cast

from apps.api.services.convex_event_relay import (
    AiidaJobEvent,
    ConvexJobProjection,
    build_convex_projection,
    compute_event_idempotency_key,
)
from apps.api.services.zpe.job_state import JobState


def make_event(**kwargs: Any) -> AiidaJobEvent:
    defaults: dict[str, Any] = {
        "job_id": "job-123",
        "node_id": "node-abc",
        "project_id": "proj-epsilon",
        "owner_id": "owner-1",
        "state": cast(JobState, "queued"),
        "event_id": "evt-1",
        "timestamp": datetime(2025, 1, 1, tzinfo=timezone.utc),
        "sequence": 1,
    }
    defaults.update(kwargs)
    return AiidaJobEvent(**defaults)


def test_build_convex_projection_maps_fields_correctly() -> None:
    event = make_event()
    projection = build_convex_projection(event)

    assert isinstance(projection, ConvexJobProjection)
    assert projection.job_id == event.job_id
    assert projection.project_id == event.project_id
    assert projection.status == event.state
    assert projection.event_time == event.timestamp
    assert projection.node_id == event.node_id
    assert projection.owner_id == event.owner_id
    assert projection.sequence == event.sequence


def test_compute_event_idempotency_key_is_stable_for_same_input() -> None:
    event = make_event()
    first_key = compute_event_idempotency_key(event)
    second_key = compute_event_idempotency_key(event)
    assert first_key == second_key


def test_compute_event_idempotency_key_changes_with_sequence() -> None:
    event = make_event()
    next_event = make_event(sequence=2)
    assert compute_event_idempotency_key(event) != compute_event_idempotency_key(next_event)


def test_projection_as_dict_matches_contract() -> None:
    event = make_event()
    projection = build_convex_projection(event)
    data = projection.as_dict()
    assert data["jobId"] == event.job_id
    assert data["projectId"] == event.project_id
    assert data["status"] == event.state
    assert data["nodeId"] == event.node_id
    assert data["ownerId"] == event.owner_id
    assert data["sequence"] == event.sequence
    assert data["eventTime"] == event.timestamp.isoformat()
