from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
import json
from typing import Any, Literal, Protocol
from urllib import error as urlerror
from urllib import request as urlrequest

from services.zpe.job_state import JobState

ProjectionStatus = Literal["queued", "running", "succeeded", "failed"]


@dataclass(frozen=True)
class AiidaJobEvent:
    """Minimal payload that arrives from AiiDA/Redis for a job projection change."""

    job_id: str
    node_id: str
    project_id: str
    owner_id: str | None
    state: JobState
    event_id: str
    timestamp: datetime
    sequence: int


@dataclass(frozen=True)
class ConvexJobProjection:
    """Field subset that GRA-16 owns in the Convex product-read projection."""

    job_id: str
    project_id: str
    status: ProjectionStatus
    event_time: datetime
    node_id: str
    owner_id: str | None
    sequence: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "jobId": self.job_id,
            "projectId": self.project_id,
            "status": self.status,
            "nodeId": self.node_id,
            "ownerId": self.owner_id,
            "sequence": self.sequence,
            "eventTime": self.event_time.isoformat(),
        }


def compute_event_idempotency_key(event: AiidaJobEvent) -> str:
    """Stable key used before persisting we already wrote this event."""

    projection_status = map_aiida_state_to_projection_status(event.state)
    digest = sha256(
        (
            f"{event.job_id}|{event.event_id}|{projection_status}|{event.sequence}"
        ).encode("utf-8")
    )
    return digest.hexdigest()


def map_aiida_state_to_projection_status(state: JobState) -> ProjectionStatus:
    state_map: dict[JobState, ProjectionStatus] = {
        "queued": "queued",
        "started": "running",
        "finished": "succeeded",
        "failed": "failed",
    }
    return state_map[state]


def build_convex_projection(event: AiidaJobEvent) -> ConvexJobProjection:
    """Project the AiiDA event onto the Convex projection schema fields."""

    return ConvexJobProjection(
        job_id=event.job_id,
        project_id=event.project_id,
        status=map_aiida_state_to_projection_status(event.state),
        event_time=event.timestamp,
        node_id=event.node_id,
        owner_id=event.owner_id,
        sequence=event.sequence,
    )


class ConvexEventDispatcher(Protocol):
    """Carries the projection payload to Convex while honoring idempotency."""

    def dispatch_job_projection(
        self,
        payload: ConvexJobProjection,
        idempotency_key: str,
    ) -> None:
        ...


class NoopConvexEventDispatcher:
    """No-op dispatcher used when relay configuration is not provided."""

    def dispatch_job_projection(
        self,
        payload: ConvexJobProjection,
        idempotency_key: str,
    ) -> None:
        _ = payload
        _ = idempotency_key


class HttpConvexEventDispatcher:
    """POST job projection updates to an HTTP relay endpoint."""

    def __init__(
        self,
        *,
        relay_url: str,
        relay_token: str | None = None,
        timeout_seconds: int = 5,
    ) -> None:
        self._relay_url = relay_url
        self._relay_token = relay_token
        self._timeout_seconds = timeout_seconds

    def dispatch_job_projection(
        self,
        payload: ConvexJobProjection,
        idempotency_key: str,
    ) -> None:
        body = json.dumps(
            {
                "projection": payload.as_dict(),
                "idempotencyKey": idempotency_key,
            },
            ensure_ascii=True,
            separators=(",", ":"),
        ).encode("utf-8")
        headers = {
            "Content-Type": "application/json",
            "Idempotency-Key": idempotency_key,
        }
        if self._relay_token:
            headers["Authorization"] = f"Bearer {self._relay_token}"
        req = urlrequest.Request(
            self._relay_url,
            method="POST",
            data=body,
            headers=headers,
        )
        try:
            with urlrequest.urlopen(req, timeout=self._timeout_seconds):
                return
        except urlerror.HTTPError as exc:
            if exc.code == 409:
                # Relay accepted an already-applied idempotent event.
                return
            raise RuntimeError(
                f"convex relay dispatch failed with HTTP {exc.code}"
            ) from exc
        except urlerror.URLError as exc:
            raise RuntimeError("convex relay dispatch failed") from exc


def get_convex_event_dispatcher(
    *,
    relay_url: str | None,
    relay_token: str | None,
    timeout_seconds: int = 5,
) -> ConvexEventDispatcher:
    if not relay_url:
        return NoopConvexEventDispatcher()
    return HttpConvexEventDispatcher(
        relay_url=relay_url,
        relay_token=relay_token,
        timeout_seconds=timeout_seconds,
    )
