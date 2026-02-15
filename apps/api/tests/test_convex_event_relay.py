from __future__ import annotations

from datetime import datetime, timezone
from email.message import Message
from io import BytesIO
from typing import Any, Literal, cast
from urllib import error as urlerror

from services.convex_event_relay import (
    AiidaJobEvent,
    ConvexJobProjection,
    HttpConvexEventDispatcher,
    build_convex_projection,
    compute_event_idempotency_key,
)
from services.zpe.job_state import JobState


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


def test_http_dispatcher_posts_projection_with_idempotency_header(
    monkeypatch: Any,
) -> None:
    event = make_event()
    projection = build_convex_projection(event)
    captured: dict[str, Any] = {}

    class _DummyResponse:
        status = 202

        def __enter__(self) -> "_DummyResponse":
            return self

        def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> Literal[False]:
            _ = (exc_type, exc, tb)
            return False

    def _fake_urlopen(request: Any, timeout: int = 0) -> _DummyResponse:
        captured["url"] = request.full_url
        captured["method"] = request.get_method()
        captured["headers"] = dict(request.header_items())
        captured["body"] = request.data.decode("utf-8")
        captured["timeout"] = timeout
        return _DummyResponse()

    monkeypatch.setattr("services.convex_event_relay.urlrequest.urlopen", _fake_urlopen)
    dispatcher = HttpConvexEventDispatcher(
        relay_url="https://relay.example.com/dispatch",
        relay_token="token-1",
        timeout_seconds=7,
    )
    dispatcher.dispatch_job_projection(projection, "idempotency-key-1")

    assert captured["url"] == "https://relay.example.com/dispatch"
    assert captured["method"] == "POST"
    assert captured["timeout"] == 7
    headers = {key.lower(): value for key, value in captured["headers"].items()}
    assert headers["idempotency-key"] == "idempotency-key-1"
    assert headers["authorization"] == "Bearer token-1"
    assert '"jobId":"job-123"' in captured["body"]
    assert '"idempotencyKey":"idempotency-key-1"' in captured["body"]


def test_http_dispatcher_treats_409_as_idempotent_success(monkeypatch: Any) -> None:
    projection = build_convex_projection(make_event())

    def _fake_urlopen(request: Any, timeout: int = 0) -> Any:
        _ = (request, timeout)
        raise urlerror.HTTPError(
            url="https://relay.example.com/dispatch",
            code=409,
            msg="conflict",
            hdrs=Message(),
            fp=BytesIO(b""),
        )

    monkeypatch.setattr("services.convex_event_relay.urlrequest.urlopen", _fake_urlopen)
    dispatcher = HttpConvexEventDispatcher(
        relay_url="https://relay.example.com/dispatch",
    )

    dispatcher.dispatch_job_projection(projection, "idempotency-key-1")
