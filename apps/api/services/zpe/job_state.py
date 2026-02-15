from __future__ import annotations

from typing import Final, Literal, TypeAlias, cast


JobState: TypeAlias = Literal["queued", "started", "finished", "failed"]

_NON_TERMINAL_STATES: Final[set[JobState]] = {"queued", "started"}
_TERMINAL_STATES: Final[set[JobState]] = {"finished", "failed"}
_ALL_STATES: Final[set[JobState]] = _NON_TERMINAL_STATES | _TERMINAL_STATES


def is_job_state(value: str) -> bool:
    return value in _ALL_STATES


def coerce_job_state(value: str) -> JobState:
    if not is_job_state(value):
        raise ValueError(f"invalid job state: {value}")
    return cast(JobState, value)


def can_transition(previous: JobState | None, next_state: JobState) -> bool:
    if previous is None:
        return True
    if previous == next_state:
        return True
    if previous == "finished":
        return False
    if previous == "failed":
        # Manual or policy-driven retry path may move failed jobs back to queue.
        return next_state == "queued"
    if previous == "queued":
        return next_state in {"started", "failed", "finished"}
    if previous == "started":
        return next_state in {"queued", "finished", "failed"}
    return False
