from __future__ import annotations

from functools import lru_cache
import json
import logging
import os
from pathlib import Path
import sqlite3
from threading import Lock
from datetime import datetime, timezone
from typing import Any
from uuid import uuid4


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def new_request_id() -> str:
    return uuid4().hex


class DurableAuditSink:
    def __init__(self, db_path: str) -> None:
        self.db_path = db_path
        self._lock = Lock()
        self._ready = False

    def _ensure_schema(self) -> None:
        if self._ready:
            return
        with self._lock:
            if self._ready:
                return
            path = Path(self.db_path)
            path.parent.mkdir(parents=True, exist_ok=True)
            with sqlite3.connect(self.db_path, timeout=5.0) as conn:
                conn.execute("PRAGMA journal_mode=WAL")
                conn.execute("PRAGMA synchronous=FULL")
                conn.execute(
                    """
                    CREATE TABLE IF NOT EXISTS zpe_audit_events (
                        event_id TEXT PRIMARY KEY,
                        occurred_at TEXT NOT NULL,
                        request_id TEXT NOT NULL,
                        tenant_id TEXT NOT NULL,
                        actor_id TEXT NOT NULL,
                        operation TEXT NOT NULL,
                        resource_type TEXT NOT NULL,
                        resource_id TEXT,
                        outcome TEXT NOT NULL,
                        payload_json TEXT NOT NULL
                    )
                    """
                )
                conn.commit()
            self._ready = True

    def write_event(
        self,
        *,
        event_id: str,
        request_id: str,
        tenant_id: str,
        actor_id: str,
        operation: str,
        resource_type: str,
        resource_id: str | None,
        outcome: str,
        metadata: dict[str, Any] | None = None,
    ) -> bool:
        fields = (
            ("event_id", event_id),
            ("request_id", request_id),
            ("tenant_id", tenant_id),
            ("actor_id", actor_id),
            ("operation", operation),
            ("resource_type", resource_type),
            ("outcome", outcome),
        )
        for name, value in fields:
            if not value or not value.strip():
                raise ValueError(f"{name} is required")
        self._ensure_schema()
        occurred_at = _now_iso()
        payload = {
            "event_id": event_id,
            "request_id": request_id,
            "tenant_id": tenant_id,
            "actor_id": actor_id,
            "operation": operation,
            "resource_type": resource_type,
            "resource_id": resource_id,
            "outcome": outcome,
            "occurred_at": occurred_at,
            "metadata": metadata or {},
        }
        payload_json = json.dumps(payload, ensure_ascii=True, separators=(",", ":"))
        with sqlite3.connect(self.db_path, timeout=5.0) as conn:
            conn.execute("PRAGMA synchronous=FULL")
            try:
                conn.execute(
                    """
                    INSERT INTO zpe_audit_events (
                        event_id,
                        occurred_at,
                        request_id,
                        tenant_id,
                        actor_id,
                        operation,
                        resource_type,
                        resource_id,
                        outcome,
                        payload_json
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        event_id,
                        occurred_at,
                        request_id,
                        tenant_id,
                        actor_id,
                        operation,
                        resource_type,
                        resource_id,
                        outcome,
                        payload_json,
                    ),
                )
                conn.commit()
                return False
            except sqlite3.IntegrityError:
                return True


@lru_cache(maxsize=1)
def get_audit_sink() -> DurableAuditSink:
    db_path = os.getenv("ZPE_AUDIT_DB_PATH", "var/zpe_audit.sqlite3")
    return DurableAuditSink(db_path=db_path)


def write_audit_event(**kwargs: Any) -> bool:
    return get_audit_sink().write_event(**kwargs)


def log_event(
    logger: logging.Logger,
    *,
    event: str,
    level: int = logging.INFO,
    **fields: Any,
) -> None:
    payload = {"ts": _now_iso(), "event": event, **fields}
    logger.log(
        level,
        json.dumps(payload, ensure_ascii=True, separators=(",", ":")),
    )
