from __future__ import annotations

import sqlite3

from services.zpe import structured_log as zpe_structured_log


def test_write_audit_event_is_idempotent(monkeypatch, tmp_path):
    db_path = tmp_path / "audit.sqlite3"
    monkeypatch.setenv("ZPE_AUDIT_DB_PATH", str(db_path))
    zpe_structured_log.get_audit_sink.cache_clear()

    first = zpe_structured_log.write_audit_event(
        event_id="event-1",
        request_id="req-1",
        tenant_id="tenant-1",
        actor_id="actor-1",
        operation="execute",
        resource_type="zpe_job",
        resource_id="job-1",
        outcome="succeeded",
        metadata={"source": "test"},
    )
    second = zpe_structured_log.write_audit_event(
        event_id="event-1",
        request_id="req-1",
        tenant_id="tenant-1",
        actor_id="actor-1",
        operation="execute",
        resource_type="zpe_job",
        resource_id="job-1",
        outcome="succeeded",
        metadata={"source": "test"},
    )
    assert first is False
    assert second is True

    with sqlite3.connect(db_path) as conn:
        row = conn.execute("SELECT COUNT(*) FROM zpe_audit_events").fetchone()
        assert row is not None
        assert row[0] == 1

