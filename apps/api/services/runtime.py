from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import logging
from pathlib import Path
import sqlite3
from threading import RLock
from typing import Any, Literal
from uuid import uuid4

from app.schemas.runtime import (
    ExecutionEvent,
    RuntimeJobStatusResponse,
    SubmitJobAccepted,
    SubmitJobCommand,
)
from app.schemas.zpe import ZPEJobStatus
from services.convex_event_relay import (
    AiidaJobEvent,
    build_convex_projection,
    compute_event_idempotency_key,
    get_convex_event_dispatcher,
)
from services.zpe.settings import get_zpe_settings

from .runtime_settings import get_runtime_settings, resolve_runtime_db_path

logger = logging.getLogger(__name__)

RuntimeState = Literal["accepted", "running", "completed", "failed"]

_STATE_ORDER: dict[RuntimeState, int] = {
    "accepted": 0,
    "running": 1,
    "completed": 2,
    "failed": 2,
}


class RuntimeConflictError(ValueError):
    pass


class RuntimeNotFoundError(KeyError):
    pass


@dataclass(frozen=True)
class SubmitResult:
    response: SubmitJobAccepted
    idempotent: bool


@dataclass(frozen=True)
class EventAck:
    idempotent: bool


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _event_time_iso(value: str) -> str:
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return _now_iso()
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.isoformat()


def _to_job_state(value: RuntimeState) -> Literal["queued", "started", "finished", "failed"]:
    if value == "accepted":
        return "queued"
    if value == "running":
        return "started"
    if value == "completed":
        return "finished"
    return "failed"


def _fingerprint_submit(command: SubmitJobCommand) -> str:
    payload = command.model_dump()
    canonical = json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
    return sha256(canonical.encode("utf-8")).hexdigest()


def _fingerprint_event(event: ExecutionEvent) -> str:
    payload = event.model_dump()
    canonical = json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
    return sha256(canonical.encode("utf-8")).hexdigest()


class RuntimeStore:
    def __init__(self, db_path: Path | None = None) -> None:
        self._db_path = db_path or resolve_runtime_db_path()
        self._lock = RLock()
        self._ensure_schema()

    def _connect(self) -> sqlite3.Connection:
        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(self._db_path), timeout=30, isolation_level=None)
        conn.row_factory = sqlite3.Row
        return conn

    def _ensure_schema(self) -> None:
        with self._connect() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS runtime_jobs (
                  job_id TEXT PRIMARY KEY,
                  tenant_id TEXT NOT NULL,
                  workspace_id TEXT NOT NULL,
                  idempotency_key TEXT NOT NULL,
                  request_fingerprint TEXT NOT NULL,
                  submission_id TEXT NOT NULL,
                  execution_owner TEXT NOT NULL,
                  state TEXT NOT NULL,
                  detail TEXT,
                  accepted_at TEXT NOT NULL,
                  updated_at TEXT NOT NULL,
                  trace_id TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE UNIQUE INDEX IF NOT EXISTS idx_runtime_submit_idempotency
                ON runtime_jobs(tenant_id, workspace_id, idempotency_key)
                """
            )
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS runtime_events (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  job_id TEXT NOT NULL,
                  event_id TEXT NOT NULL,
                  event_fingerprint TEXT NOT NULL,
                  occurred_at TEXT NOT NULL,
                  state TEXT NOT NULL,
                  detail TEXT,
                  trace_id TEXT NOT NULL,
                  created_at TEXT NOT NULL,
                  UNIQUE(job_id, event_id)
                )
                """
            )

    def submit(self, command: SubmitJobCommand, *, trace_id: str) -> SubmitResult:
        fingerprint = _fingerprint_submit(command)
        now = _now_iso()
        with self._lock, self._connect() as conn:
            row = conn.execute(
                """
                SELECT * FROM runtime_jobs
                WHERE tenant_id = ? AND workspace_id = ? AND idempotency_key = ?
                """,
                (command.tenant_id, command.workspace_id, command.idempotency_key),
            ).fetchone()
            if row is not None:
                if row["request_fingerprint"] != fingerprint:
                    raise RuntimeConflictError("idempotency_conflict")
                return SubmitResult(
                    response=SubmitJobAccepted(
                        job_id=row["job_id"],
                        submission_id=row["submission_id"],
                        execution_owner=row["execution_owner"],
                        accepted_at=row["accepted_at"],
                        trace_id=row["trace_id"],
                    ),
                    idempotent=True,
                )

            existing_job = conn.execute(
                "SELECT job_id FROM runtime_jobs WHERE job_id = ?",
                (command.job_id,),
            ).fetchone()
            if existing_job is not None:
                raise RuntimeConflictError("job_id_conflict")

            submission_id = f"sub-{uuid4().hex}"
            execution_owner = get_runtime_settings().execution_owner
            conn.execute(
                """
                INSERT INTO runtime_jobs(
                  job_id, tenant_id, workspace_id, idempotency_key, request_fingerprint,
                  submission_id, execution_owner, state, detail, accepted_at, updated_at, trace_id
                ) VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    command.job_id,
                    command.tenant_id,
                    command.workspace_id,
                    command.idempotency_key,
                    fingerprint,
                    submission_id,
                    execution_owner,
                    "accepted",
                    None,
                    now,
                    now,
                    trace_id,
                ),
            )
            return SubmitResult(
                response=SubmitJobAccepted(
                    job_id=command.job_id,
                    submission_id=submission_id,
                    execution_owner=execution_owner,
                    accepted_at=now,
                    trace_id=trace_id,
                ),
                idempotent=False,
            )

    def apply_event(self, event: ExecutionEvent) -> EventAck:
        fingerprint = _fingerprint_event(event)
        now = _now_iso()
        with self._lock, self._connect() as conn:
            job = conn.execute(
                "SELECT * FROM runtime_jobs WHERE job_id = ?",
                (event.job_id,),
            ).fetchone()
            if job is None:
                raise RuntimeNotFoundError(event.job_id)

            if (
                job["tenant_id"] != event.tenant_id
                or job["workspace_id"] != event.workspace_id
                or job["submission_id"] != event.submission_id
            ):
                raise RuntimeConflictError("scope_conflict")

            existing_event = conn.execute(
                "SELECT event_fingerprint FROM runtime_events WHERE job_id = ? AND event_id = ?",
                (event.job_id, event.event_id),
            ).fetchone()
            if existing_event is not None:
                if existing_event["event_fingerprint"] != fingerprint:
                    raise RuntimeConflictError("event_id_conflict")
                return EventAck(idempotent=True)

            prev_state = job["state"]
            next_state = event.state
            if _STATE_ORDER[next_state] < _STATE_ORDER[prev_state]:
                raise RuntimeConflictError("invalid_transition")
            if prev_state in {"completed", "failed"} and next_state != prev_state:
                raise RuntimeConflictError("invalid_transition")

            conn.execute(
                """
                INSERT INTO runtime_events(
                  job_id, event_id, event_fingerprint, occurred_at, state, detail, trace_id, created_at
                ) VALUES(?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    event.job_id,
                    event.event_id,
                    fingerprint,
                    _event_time_iso(event.occurred_at),
                    event.state,
                    event.status_detail,
                    event.trace_id,
                    now,
                ),
            )
            conn.execute(
                """
                UPDATE runtime_jobs
                SET state = ?, detail = ?, updated_at = ?, trace_id = ?
                WHERE job_id = ?
                """,
                (
                    next_state,
                    event.status_detail,
                    _event_time_iso(event.occurred_at),
                    event.trace_id,
                    event.job_id,
                ),
            )

        self._dispatch_projection(event)
        return EventAck(idempotent=False)

    def _dispatch_projection(self, event: ExecutionEvent) -> None:
        settings = get_zpe_settings()
        dispatcher = get_convex_event_dispatcher(
            relay_url=settings.convex_relay_url,
            relay_token=settings.convex_relay_token,
            timeout_seconds=settings.convex_relay_timeout_seconds,
        )
        try:
            evt = AiidaJobEvent(
                job_id=event.job_id,
                node_id=event.execution_id,
                project_id=event.workspace_id,
                owner_id=None,
                state=_to_job_state(event.state),
                event_id=event.event_id,
                timestamp=datetime.fromisoformat(_event_time_iso(event.occurred_at)),
                sequence=1,
            )
            projection = build_convex_projection(evt)
            dispatcher.dispatch_job_projection(
                payload=projection,
                idempotency_key=compute_event_idempotency_key(evt),
            )
        except Exception:
            logger.warning("runtime projection dispatch failed", exc_info=True)

    def get_status(self, job_id: str, *, tenant_id: str) -> RuntimeJobStatusResponse:
        with self._connect() as conn:
            row = conn.execute(
                "SELECT * FROM runtime_jobs WHERE job_id = ? AND tenant_id = ?",
                (job_id, tenant_id),
            ).fetchone()
        if row is None:
            raise RuntimeNotFoundError(job_id)
        return RuntimeJobStatusResponse(
            job_id=row["job_id"],
            submission_id=row["submission_id"],
            execution_owner=row["execution_owner"],
            state=row["state"],
            detail=row["detail"],
            updated_at=row["updated_at"],
            trace_id=row["trace_id"],
        )

    def get_projection_status(self, job_id: str, *, tenant_id: str) -> ZPEJobStatus:
        status = self.get_status(job_id, tenant_id=tenant_id)
        return ZPEJobStatus(
            status=_to_job_state(status.state),
            detail=status.detail,
            updated_at=status.updated_at,
        )

    def get_detail(self, job_id: str, *, tenant_id: str) -> dict[str, Any]:
        with self._connect() as conn:
            job = conn.execute(
                "SELECT * FROM runtime_jobs WHERE job_id = ? AND tenant_id = ?",
                (job_id, tenant_id),
            ).fetchone()
            event = conn.execute(
                """
                SELECT event_id, occurred_at, state, detail, trace_id
                FROM runtime_events
                WHERE job_id = ?
                ORDER BY id DESC
                LIMIT 1
                """,
                (job_id,),
            ).fetchone()
        if job is None:
            raise RuntimeNotFoundError(job_id)
        return {
            "job_id": job["job_id"],
            "submission_id": job["submission_id"],
            "state": job["state"],
            "detail": job["detail"],
            "updated_at": job["updated_at"],
            "last_event": dict(event) if event is not None else None,
        }


_store_instance: RuntimeStore | None = None


def get_runtime_store() -> RuntimeStore:
    global _store_instance
    if _store_instance is None:
        _store_instance = RuntimeStore()
    return _store_instance
