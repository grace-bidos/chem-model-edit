from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from hashlib import sha256
from typing import Any, Protocol

from services.zpe.job_state import JobState


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
    status: JobState
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

    digest = sha256(
        f"{event.job_id}|{event.event_id}|{event.state}|{event.sequence}".encode("utf-8")
    )
    return digest.hexdigest()


def build_convex_projection(event: AiidaJobEvent) -> ConvexJobProjection:
    """Project the AiiDA event onto the Convex projection schema fields."""

    return ConvexJobProjection(
        job_id=event.job_id,
        project_id=event.project_id,
        status=event.state,
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
