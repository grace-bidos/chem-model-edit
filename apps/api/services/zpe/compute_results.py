from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import logging
from typing import Any, Dict, Optional, cast

from redis import Redis, WatchError

from services.convex_event_relay import (
    AiidaJobEvent,
    build_convex_projection,
    compute_event_idempotency_key,
    get_convex_event_dispatcher,
)

from .job_meta import get_job_meta_store
from .job_owner import get_job_owner_store
from .job_state import JobState, can_transition, coerce_job_state
from .queue import get_redis_connection
from .settings import get_zpe_settings


_STATUS_PREFIX = "zpe:status:"
_RESULT_PREFIX = "zpe:result:"
_SUMMARY_PREFIX = "zpe:summary:"
_FREQS_PREFIX = "zpe:freqs:"
_LEASE_PREFIX = "zpe:lease:"
_LEASE_INDEX = "zpe:lease:index"
_RETRY_PREFIX = "zpe:retry_count:"
_DELAY_ZSET = "zpe:delay"
_DLQ = "zpe:dlq"
_FAILURE_SUBMIT_PREFIX = "zpe:failure_submit:"
_RELAY_DISPATCH_PREFIX = "zpe:convex:relay:dispatch:"
_RELAY_SEQUENCE_FIELD = "relay_sequence"


logger = logging.getLogger(__name__)


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _decode_map(data: Dict[bytes, bytes]) -> Dict[str, str]:
    return {k.decode("utf-8"): v.decode("utf-8") for k, v in data.items()}


def _get_lease(redis: Redis, job_id: str) -> Dict[str, str]:
    raw = cast(dict[bytes, bytes], redis.hgetall(f"{_LEASE_PREFIX}{job_id}"))
    return _decode_map(raw) if raw else {}


def _decode_status(raw: Optional[str]) -> JobState | None:
    if raw is None:
        return None
    return coerce_job_state(raw)


def _to_iso_datetime(value: str | None) -> datetime:
    if not value:
        return datetime.now(timezone.utc)
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return datetime.now(timezone.utc)
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed


def _coerce_non_empty(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        return stripped or None
    rendered = str(value).strip()
    return rendered or None


def _dispatch_runtime_state_transition(
    *,
    redis: Redis,
    job_id: str,
    state: JobState,
    sequence: int,
    updated_at: str,
    ttl_seconds: int,
) -> None:
    if sequence <= 0:
        return
    dispatch_key = f"{_RELAY_DISPATCH_PREFIX}{job_id}:{sequence}"
    acquired = cast(
        bool,
        redis.set(
            dispatch_key,
            "1",
            nx=True,
            ex=ttl_seconds,
        ),
    )
    if not acquired:
        return
    settings = get_zpe_settings()
    dispatcher = get_convex_event_dispatcher(
        relay_url=settings.convex_relay_url,
        relay_token=settings.convex_relay_token,
        timeout_seconds=settings.convex_relay_timeout_seconds,
    )
    owner_id = get_job_owner_store().get_owner(job_id)
    meta = get_job_meta_store().get_meta(job_id)
    project_id = _coerce_non_empty(meta.get("project_id")) or _coerce_non_empty(
        meta.get("queue_name")
    )
    node_id = _coerce_non_empty(meta.get("node_id")) or job_id
    event = AiidaJobEvent(
        job_id=job_id,
        node_id=node_id,
        project_id=project_id or "zpe",
        owner_id=owner_id,
        state=state,
        event_id=f"runtime:{job_id}:{sequence}",
        timestamp=_to_iso_datetime(updated_at),
        sequence=sequence,
    )
    projection = build_convex_projection(event)
    idempotency_key = compute_event_idempotency_key(event)
    try:
        dispatcher.dispatch_job_projection(
            payload=projection,
            idempotency_key=idempotency_key,
        )
    except Exception:
        redis.delete(dispatch_key)
        logger.warning(
            "convex relay dispatch failed for runtime transition",
            extra={
                "job_id": job_id,
                "state": state,
                "sequence": sequence,
            },
            exc_info=True,
        )


@dataclass
class ResultSubmitOutcome:
    idempotent: bool


@dataclass
class FailureOutcome:
    requeued: bool
    retry_count: int


def _validate_recorded_failure_payload(
    *,
    existing: Dict[str, str],
    worker_id: str,
    error_code: str,
    error_message: str,
    traceback: Optional[str],
) -> None:
    if existing.get("worker_id") != worker_id:
        raise PermissionError("lease mismatch")
    if (
        existing.get("error_code") != error_code
        or existing.get("error_message") != error_message
        or existing.get("traceback", "") != (traceback or "")
    ):
        raise ValueError("failure already submitted with different payload")


def _failure_outcome_from_record(existing: Dict[str, str]) -> FailureOutcome:
    return FailureOutcome(
        requeued=existing.get("requeued") == "1",
        retry_count=int(existing.get("retry_count", "0")),
    )


def submit_result(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    result: Dict[str, Any],
    summary_text: str,
    freqs_csv: str,
) -> ResultSubmitOutcome:
    settings = get_zpe_settings()
    redis = get_redis_connection()
    status_key = f"{_STATUS_PREFIX}{job_id}"
    result_key = f"{_RESULT_PREFIX}{job_id}"
    summary_key = f"{_SUMMARY_PREFIX}{job_id}"
    freqs_key = f"{_FREQS_PREFIX}{job_id}"
    lease_key = f"{_LEASE_PREFIX}{job_id}"
    ttl = settings.result_ttl_seconds
    payload_json = json.dumps(result)

    for _ in range(3):
        pipe = redis.pipeline()
        pipe_any = cast(Any, pipe)
        try:
            pipe_any.watch(status_key, lease_key)
            status = _decode_map(cast(dict[bytes, bytes], pipe.hgetall(status_key)))
            previous = _decode_status(status.get("status"))
            if previous == "finished":
                existing_result = pipe.get(result_key)
                existing_summary = pipe.get(summary_key)
                existing_freqs = pipe.get(freqs_key)
                pipe.unwatch()
                if (
                    existing_result == payload_json.encode("utf-8")
                    and existing_summary == summary_text.encode("utf-8")
                    and existing_freqs == freqs_csv.encode("utf-8")
                ):
                    sequence = int(status.get(_RELAY_SEQUENCE_FIELD) or "1")
                    _dispatch_runtime_state_transition(
                        redis=redis,
                        job_id=job_id,
                        state="finished",
                        sequence=sequence,
                        updated_at=status.get("updated_at") or _now_iso(),
                        ttl_seconds=ttl,
                    )
                    return ResultSubmitOutcome(idempotent=True)
                raise ValueError("result already submitted")
            if not can_transition(previous, "finished"):
                raise ValueError(
                    f"invalid job state transition: {previous or '<none>'} -> finished"
                )

            lease = _get_lease(pipe, job_id)
            if not lease:
                raise PermissionError("lease not found")
            if lease.get("worker_id") != worker_id or lease.get("lease_id") != lease_id:
                raise PermissionError("lease mismatch")

            updated_at = _now_iso()
            pipe.multi()
            pipe.setex(result_key, ttl, payload_json)
            pipe.setex(summary_key, ttl, summary_text)
            pipe.setex(freqs_key, ttl, freqs_csv)
            pipe.hset(
                status_key,
                mapping={"status": "finished", "detail": "", "updated_at": updated_at},
            )
            pipe.hincrby(status_key, _RELAY_SEQUENCE_FIELD, 1)
            pipe.expire(status_key, ttl)
            pipe.delete(lease_key)
            pipe.zrem(_LEASE_INDEX, job_id)
            response = cast(list[Any], pipe.execute())
            sequence = int(response[4])
            _dispatch_runtime_state_transition(
                redis=redis,
                job_id=job_id,
                state="finished",
                sequence=sequence,
                updated_at=updated_at,
                ttl_seconds=ttl,
            )
            return ResultSubmitOutcome(idempotent=False)
        except WatchError:
            continue
        finally:
            pipe.reset()
    raise RuntimeError("failed to submit result due to concurrent updates")


def submit_failure(
    *,
    job_id: str,
    worker_id: str,
    lease_id: str,
    error_code: str,
    error_message: str,
    traceback: Optional[str] = None,
) -> FailureOutcome:
    settings = get_zpe_settings()
    redis = get_redis_connection()
    status_key = f"{_STATUS_PREFIX}{job_id}"
    lease_key = f"{_LEASE_PREFIX}{job_id}"
    retry_key = f"{_RETRY_PREFIX}{job_id}"
    submit_key = f"{_FAILURE_SUBMIT_PREFIX}{job_id}:{lease_id}"
    ttl = settings.result_ttl_seconds

    existing = _decode_map(cast(dict[bytes, bytes], redis.hgetall(submit_key)))
    if existing:
        _validate_recorded_failure_payload(
            existing=existing,
            worker_id=worker_id,
            error_code=error_code,
            error_message=error_message,
            traceback=traceback,
        )
        lease = _get_lease(redis, job_id)
        if not lease:
            return _failure_outcome_from_record(existing)

    for _ in range(3):
        pipe = redis.pipeline()
        pipe_any = cast(Any, pipe)
        try:
            pipe_any.watch(status_key, lease_key, retry_key, submit_key)
            existing = _decode_map(cast(dict[bytes, bytes], pipe.hgetall(submit_key)))
            if existing:
                _validate_recorded_failure_payload(
                    existing=existing,
                    worker_id=worker_id,
                    error_code=error_code,
                    error_message=error_message,
                    traceback=traceback,
                )
            status = _decode_map(cast(dict[bytes, bytes], pipe.hgetall(status_key)))
            previous = _decode_status(status.get("status"))
            lease = _get_lease(pipe, job_id)
            if existing and not lease:
                pipe.unwatch()
                return _failure_outcome_from_record(existing)
            if not lease:
                existing = _decode_map(cast(dict[bytes, bytes], pipe.hgetall(submit_key)))
                if existing:
                    _validate_recorded_failure_payload(
                        existing=existing,
                        worker_id=worker_id,
                        error_code=error_code,
                        error_message=error_message,
                        traceback=traceback,
                    )
                    pipe.unwatch()
                    return _failure_outcome_from_record(existing)
                raise PermissionError("lease not found")
            if lease.get("worker_id") != worker_id or lease.get("lease_id") != lease_id:
                raise PermissionError("lease mismatch")

            current_retry = cast(Optional[bytes], pipe.get(retry_key))
            retry_count = int(current_retry or 0) + 1
            detail = f"{error_code}: {error_message}"
            next_state: JobState = (
                "queued" if retry_count <= settings.retry_max else "failed"
            )
            if not can_transition(previous, next_state):
                raise ValueError(
                    f"invalid job state transition: {previous or '<none>'} -> {next_state}"
                )

            updated_at = _now_iso()
            pipe.multi()
            pipe.set(retry_key, retry_count)
            if retry_count <= settings.retry_max:
                delay = min(
                    settings.retry_max_delay_seconds,
                    settings.retry_base_delay_seconds * (2 ** (retry_count - 1)),
                )
                run_at = int(datetime.now(timezone.utc).timestamp()) + int(delay)
                pipe.zadd(_DELAY_ZSET, {job_id: run_at})
                pipe.hset(
                    status_key,
                    mapping={
                        "status": "queued",
                        "detail": f"requeue in {delay}s ({detail})",
                        "updated_at": updated_at,
                    },
                )
                pipe.hincrby(status_key, _RELAY_SEQUENCE_FIELD, 1)
                pipe.expire(status_key, ttl)
                pipe.expire(retry_key, ttl)
                pipe.delete(lease_key)
                pipe.zrem(_LEASE_INDEX, job_id)
                pipe.hset(
                    submit_key,
                    mapping={
                        "worker_id": worker_id,
                        "error_code": error_code,
                        "error_message": error_message,
                        "traceback": traceback or "",
                        "retry_count": str(retry_count),
                        "requeued": "1",
                    },
                )
                pipe.expire(submit_key, ttl)
                response = cast(list[Any], pipe.execute())
                sequence = int(response[3])
                _dispatch_runtime_state_transition(
                    redis=redis,
                    job_id=job_id,
                    state="queued",
                    sequence=sequence,
                    updated_at=updated_at,
                    ttl_seconds=ttl,
                )
                return FailureOutcome(requeued=True, retry_count=retry_count)

            pipe.lpush(_DLQ, job_id)
            pipe.hset(
                status_key,
                mapping={
                    "status": "failed",
                    "detail": detail,
                    "updated_at": updated_at,
                },
            )
            pipe.hincrby(status_key, _RELAY_SEQUENCE_FIELD, 1)
            pipe.expire(status_key, ttl)
            pipe.expire(retry_key, ttl)
            pipe.delete(lease_key)
            pipe.zrem(_LEASE_INDEX, job_id)
            pipe.hset(
                submit_key,
                mapping={
                    "worker_id": worker_id,
                    "error_code": error_code,
                    "error_message": error_message,
                    "traceback": traceback or "",
                    "retry_count": str(retry_count),
                    "requeued": "0",
                },
            )
            pipe.expire(submit_key, ttl)
            response = cast(list[Any], pipe.execute())
            sequence = int(response[3])
            _dispatch_runtime_state_transition(
                redis=redis,
                job_id=job_id,
                state="failed",
                sequence=sequence,
                updated_at=updated_at,
                ttl_seconds=ttl,
            )
            return FailureOutcome(requeued=False, retry_count=retry_count)
        except WatchError:
            continue
        finally:
            pipe.reset()
    raise RuntimeError("failed to submit failure due to concurrent updates")
