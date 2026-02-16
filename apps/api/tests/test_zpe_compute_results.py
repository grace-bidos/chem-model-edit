from __future__ import annotations

import fakeredis
import pytest

from services.convex_event_relay import ConvexJobProjection
from services.zpe import compute_results as zpe_results
from services.zpe.settings import ZPESettings


class _StaticOwnerStore:
    def __init__(self, owner_id: str | None) -> None:
        self._owner_id = owner_id

    def get_owner(self, job_id: str) -> str | None:
        _ = job_id
        return self._owner_id


class _StaticMetaStore:
    def __init__(self, meta: dict[str, object]) -> None:
        self._meta = meta

    def get_meta(self, job_id: str) -> dict[str, object]:
        _ = job_id
        return self._meta


class _RecordingDispatcher:
    def __init__(self, *, fail_first: bool = False) -> None:
        self.fail_first = fail_first
        self.calls: list[tuple[ConvexJobProjection, str]] = []

    def dispatch_job_projection(
        self,
        payload: ConvexJobProjection,
        idempotency_key: str,
    ) -> None:
        self.calls.append((payload, idempotency_key))
        if self.fail_first and len(self.calls) == 1:
            raise RuntimeError("transient relay failure")


def _patch(monkeypatch, *, dispatcher: _RecordingDispatcher | None = None):
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(zpe_results, "get_redis_connection", lambda: fake)
    settings = ZPESettings(
        result_ttl_seconds=60,
        retry_max=2,
        retry_base_delay_seconds=1,
        retry_max_delay_seconds=5,
        convex_relay_url="https://relay.example.com/dispatch",
    )
    monkeypatch.setattr(zpe_results, "get_zpe_settings", lambda: settings)
    used_dispatcher = dispatcher or _RecordingDispatcher()
    monkeypatch.setattr(
        zpe_results,
        "get_convex_event_dispatcher",
        lambda **kwargs: used_dispatcher,
    )
    monkeypatch.setattr(
        zpe_results,
        "get_job_owner_store",
        lambda: _StaticOwnerStore("owner-1"),
    )
    monkeypatch.setattr(
        zpe_results,
        "get_job_meta_store",
        lambda: _StaticMetaStore({"queue_name": "zpe"}),
    )
    return fake, used_dispatcher


def test_submit_result_idempotent(monkeypatch):
    fake, dispatcher = _patch(monkeypatch)
    job_id = "job-1"
    lease_id = "lease-1"
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    outcome = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.1},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert outcome.idempotent is False

    outcome_again = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.1},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert outcome_again.idempotent is True
    assert len(dispatcher.calls) == 1
    assert dispatcher.calls[0][0].status == "finished"
    assert dispatcher.calls[0][0].sequence == 1


def test_submit_result_duplicate_replay_retries_dispatch(monkeypatch):
    dispatcher = _RecordingDispatcher(fail_first=True)
    fake, _ = _patch(monkeypatch, dispatcher=dispatcher)
    job_id = "job-1b"
    lease_id = "lease-1b"
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    first = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.2},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert first.idempotent is False

    replay = zpe_results.submit_result(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        result={"zpe_ev": 0.2},
        summary_text="summary",
        freqs_csv="freqs",
    )
    assert replay.idempotent is True
    assert len(dispatcher.calls) == 2


def test_submit_failure_requeue_and_dlq(monkeypatch):
    fake, dispatcher = _patch(monkeypatch)
    job_id = "job-2"
    lease_id = "lease-2"
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    outcome = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome.requeued is True
    assert outcome.retry_count == 1
    assert fake.zcard("zpe:delay") == 1

    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )
    outcome2 = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome2.requeued is True
    assert outcome2.retry_count == 2

    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )
    outcome3 = zpe_results.submit_failure(
        job_id=job_id,
        worker_id="worker-1",
        lease_id=lease_id,
        error_code="ERR",
        error_message="boom",
    )
    assert outcome3.requeued is False
    assert outcome3.retry_count == 3
    assert fake.lrange("zpe:dlq", 0, -1) == [job_id.encode("utf-8")]
    assert [call[0].status for call in dispatcher.calls] == ["queued", "queued", "failed"]
    assert [call[0].sequence for call in dispatcher.calls] == [1, 2, 3]


def test_submit_failure_rejects_transition_from_finished(monkeypatch):
    fake, _dispatcher = _patch(monkeypatch)
    job_id = "job-finished"
    lease_id = "lease-finished"
    fake.hset(
        f"zpe:status:{job_id}",
        mapping={"status": "finished", "detail": "", "updated_at": "2026-01-01T00:00:00+00:00"},
    )
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    with pytest.raises(ValueError, match="invalid job state transition"):
        zpe_results.submit_failure(
            job_id=job_id,
            worker_id="worker-1",
            lease_id=lease_id,
            error_code="ERR",
            error_message="boom",
        )


def test_submit_result_rejects_transition_from_failed(monkeypatch):
    fake, _dispatcher = _patch(monkeypatch)
    job_id = "job-failed"
    lease_id = "lease-failed"
    fake.hset(
        f"zpe:status:{job_id}",
        mapping={
            "status": "failed",
            "detail": "boom",
            "updated_at": "2026-01-01T00:00:00+00:00",
        },
    )
    fake.hset(
        f"zpe:lease:{job_id}",
        mapping={"worker_id": "worker-1", "lease_id": lease_id, "expires_at": "x"},
    )

    with pytest.raises(ValueError, match="invalid job state transition"):
        zpe_results.submit_result(
            job_id=job_id,
            worker_id="worker-1",
            lease_id=lease_id,
            result={"zpe_ev": 0.1},
            summary_text="summary",
            freqs_csv="freqs",
        )
