from __future__ import annotations

import json
import logging
import socket
import time
import traceback
from datetime import datetime, timezone
from typing import Any, Dict, Optional, Tuple
from urllib import request as urlrequest
from urllib import error as urlerror

from .settings import get_zpe_settings
from .structured_log import log_event
from .worker import compute_zpe_artifacts

logger = logging.getLogger("zpe.http_worker")


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _request_json(
    method: str,
    url: str,
    *,
    token: str,
    payload: Optional[Dict[str, Any]] = None,
    timeout: int = 60,
) -> Tuple[int, Optional[Dict[str, Any]]]:
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }
    data = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(payload).encode("utf-8")
    req = urlrequest.Request(url, data=data, headers=headers, method=method)
    try:
        with urlrequest.urlopen(req, timeout=timeout) as resp:
            status = resp.status
            body = resp.read().decode("utf-8")
    except urlerror.HTTPError as exc:
        status = exc.code
        body = exc.read().decode("utf-8") if exc.fp else ""
    except urlerror.URLError:
        return 0, None
    if not body:
        return status, None
    try:
        return status, json.loads(body)
    except json.JSONDecodeError:
        return status, None


def _lease_job(base_url: str, token: str, timeout: int) -> Optional[Dict[str, Any]]:
    status, data = _request_json(
        "POST",
        f"{base_url}/api/zpe/compute/jobs/lease",
        token=token,
        payload=None,
        timeout=timeout,
    )
    if status == 204:
        return None
    if status != 200:
        raise RuntimeError(f"lease failed (status={status})")
    return data


def _submit_result(
    base_url: str,
    token: str,
    job_id: str,
    tenant_id: str,
    lease_id: str,
    result: Dict[str, Any],
    summary_text: str,
    freqs_csv: str,
    meta: Dict[str, Any],
    execution_event: Dict[str, Any] | None = None,
    timeout: int = 60,
) -> None:
    payload = {
        "tenant_id": tenant_id,
        "lease_id": lease_id,
        "result": result,
        "summary_text": summary_text,
        "freqs_csv": freqs_csv,
        "meta": meta,
    }
    if execution_event is not None:
        payload["execution_event"] = execution_event
    status, _ = _request_json(
        "POST",
        f"{base_url}/api/zpe/compute/jobs/{job_id}/result",
        token=token,
        payload=payload,
        timeout=timeout,
    )
    if status not in {200, 409}:
        raise RuntimeError(f"result submit failed (status={status})")


def _submit_failed(
    base_url: str,
    token: str,
    job_id: str,
    tenant_id: str,
    lease_id: str,
    error_code: str,
    error_message: str,
    tb: str,
    execution_event: Dict[str, Any] | None = None,
    timeout: int = 60,
) -> None:
    payload = {
        "tenant_id": tenant_id,
        "lease_id": lease_id,
        "error_code": error_code,
        "error_message": error_message,
        "traceback": tb,
    }
    if execution_event is not None:
        payload["execution_event"] = execution_event
    _request_json(
        "POST",
        f"{base_url}/api/zpe/compute/jobs/{job_id}/failed",
        token=token,
        payload=payload,
        timeout=timeout,
    )


def run_http_worker() -> None:
    settings = get_zpe_settings()
    if not settings.worker_token:
        raise RuntimeError("ZPE_WORKER_TOKEN is required for HTTP worker")

    base_url = settings.control_api_url.rstrip("/")
    if base_url.endswith("/api"):
        base_url = base_url[: -len("/api")]
    token = settings.worker_token
    poll_interval = max(1, int(settings.worker_poll_interval_seconds))
    max_interval = max(poll_interval, int(settings.worker_poll_max_interval_seconds))
    timeout = max(1, int(settings.worker_request_timeout_seconds))

    backoff = poll_interval
    hostname = socket.gethostname()
    log_event(
        logger,
        event="zpe_worker_start",
        service="compute-plane",
        stage="worker",
        status="started",
        worker_hostname=hostname,
        backend="remote-http",
        result_store=settings.result_store,
    )

    while True:
        try:
            lease = _lease_job(base_url, token, timeout)
        except Exception as exc:
            log_event(
                logger,
                event="zpe_lease_error",
                service="compute-plane",
                stage="lease",
                status="error",
                worker_hostname=hostname,
                backend="remote-http",
                result_store=settings.result_store,
                error_message=str(exc),
            )
            time.sleep(backoff)
            backoff = min(max_interval, backoff * 2)
            continue

        if lease is None:
            time.sleep(backoff)
            backoff = min(max_interval, backoff * 2)
            continue

        backoff = poll_interval
        job_id = lease.get("job_id")
        payload = lease.get("payload")
        lease_id = lease.get("lease_id")
        meta = lease.get("meta") or {}
        tenant_id = meta.get("tenant_id")
        request_id = meta.get("request_id")
        user_id = meta.get("user_id")
        if not job_id or not payload or not lease_id or not isinstance(tenant_id, str):
            log_event(
                logger,
                event="zpe_lease_invalid",
                service="compute-plane",
                stage="lease",
                status="invalid",
                worker_hostname=hostname,
                backend="remote-http",
                result_store=settings.result_store,
            )
            time.sleep(poll_interval)
            continue

        start = time.monotonic()
        trace_id = str(request_id or f"trace-{job_id}-{lease_id}")
        workspace_id = str(meta.get("workspace_id") or tenant_id)
        submission_id = str(meta.get("submission_id") or request_id or lease_id)
        management_node_id = str(meta.get("management_node_id") or hostname)
        execution_id = str(
            meta.get("execution_id") or f"{management_node_id}:{job_id}:{lease_id}"
        )
        scheduler_ref: Dict[str, Any] = {}
        for key, source in (
            ("slurm_job_id", "slurm_job_id"),
            ("partition", "slurm_partition"),
            ("qos", "slurm_qos"),
        ):
            value = meta.get(source)
            if isinstance(value, str) and value.strip():
                scheduler_ref[key] = value.strip()
        try:
            log_event(
                logger,
                event="zpe_compute_start",
                service="compute-plane",
                stage="compute",
                status="started",
                request_id=request_id,
                job_id=job_id,
                user_id=user_id,
                worker_hostname=hostname,
                backend="remote-http",
                result_store=settings.result_store,
            )
            artifacts = compute_zpe_artifacts(payload, job_id=job_id)
            elapsed = time.monotonic() - start
            meta = {
                "worker_hostname": hostname,
                "computation_time_seconds": elapsed,
                "timestamp": _now_iso(),
                "request_id": request_id,
            }
            _submit_result(
                base_url,
                token,
                job_id,
                tenant_id,
                lease_id,
                artifacts.result,
                artifacts.summary_text,
                artifacts.freqs_csv,
                meta,
                execution_event={
                    "event_id": f"execution:{job_id}:{lease_id}:completed",
                    "tenant_id": tenant_id,
                    "workspace_id": workspace_id,
                    "job_id": job_id,
                    "submission_id": submission_id,
                    "execution_id": execution_id,
                    "state": "completed",
                    "occurred_at": _now_iso(),
                    "trace_id": trace_id,
                    "status_detail": "job completed",
                    "scheduler_ref": scheduler_ref or None,
                    "result_ref": {
                        "output_uri": f"zpe://jobs/{job_id}/result",
                        "metadata_uri": f"zpe://jobs/{job_id}/files/summary",
                    },
                },
                timeout=timeout,
            )
            log_event(
                logger,
                event="zpe_compute_success",
                service="compute-plane",
                stage="compute",
                status="success",
                request_id=request_id,
                job_id=job_id,
                user_id=user_id,
                worker_hostname=hostname,
                backend="remote-http",
                result_store=settings.result_store,
                duration_ms=int(elapsed * 1000),
                exit_code=0,
                qe_version=None,
            )
        except Exception as exc:
            tb = traceback.format_exc()
            elapsed = time.monotonic() - start
            log_event(
                logger,
                event="zpe_compute_failed",
                service="compute-plane",
                stage="compute",
                status="failed",
                request_id=request_id,
                job_id=job_id,
                user_id=user_id,
                worker_hostname=hostname,
                backend="remote-http",
                result_store=settings.result_store,
                duration_ms=int(elapsed * 1000),
                exit_code=1,
                qe_version=None,
                error_message=str(exc),
            )
            _submit_failed(
                base_url,
                token,
                job_id,
                tenant_id,
                lease_id,
                "COMPUTE_ERROR",
                str(exc),
                tb,
                execution_event={
                    "event_id": f"execution:{job_id}:{lease_id}:failed",
                    "tenant_id": tenant_id,
                    "workspace_id": workspace_id,
                    "job_id": job_id,
                    "submission_id": submission_id,
                    "execution_id": execution_id,
                    "state": "failed",
                    "occurred_at": _now_iso(),
                    "trace_id": trace_id,
                    "status_detail": str(exc),
                    "scheduler_ref": scheduler_ref or None,
                    "error": {
                        "code": "COMPUTE_ERROR",
                        "message": str(exc),
                        "retryable": True,
                    },
                },
                timeout=timeout,
            )


if __name__ == "__main__":
    run_http_worker()
