from __future__ import annotations

import json
from http.server import BaseHTTPRequestHandler, HTTPServer
import threading
from typing import Any, Dict

from app.schemas.zpe import ComputeFailedRequest, ComputeResultRequest
from services.zpe import http_worker
from services.zpe import worker as zpe_worker
from services.zpe.settings import ZPESettings

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=1, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
""".strip()


def _make_handler() -> type[BaseHTTPRequestHandler]:
    class _LeaseHandler(BaseHTTPRequestHandler):
        token = "test-token"  # noqa: S105
        job_id = "job-1"
        lease_id = "lease-1"
        leased = False
        result_payload: Dict[str, Any] | None = None
        failed_payload: Dict[str, Any] | None = None

        def _send_json(
            self, status: int, payload: Dict[str, Any] | None = None
        ) -> None:
            self.send_response(status)
            if payload is not None:
                body = json.dumps(payload).encode("utf-8")
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
            else:
                self.end_headers()

        def _auth_ok(self) -> bool:
            auth = self.headers.get("Authorization")
            return auth == f"Bearer {self.token}"

        def do_POST(self) -> None:
            handler_cls = type(self)
            if not self._auth_ok():
                self._send_json(401, {"detail": "unauthorized"})
                return

            if self.path == "/api/zpe/compute/jobs/lease":
                if handler_cls.leased:
                    self._send_json(204)
                    return
                handler_cls.leased = True
                payload = {
                    "job_id": self.job_id,
                    "lease_id": self.lease_id,
                    "lease_ttl_seconds": 30,
                    "meta": {
                        "tenant_id": "tenant-1",
                        "request_id": "trace-req-1",
                        "workspace_id": "workspace-1",
                        "submission_id": "submission-1",
                        "management_node_id": "node-1",
                        "slurm_job_id": "12345",
                        "slurm_partition": "short",
                        "slurm_qos": "normal",
                    },
                    "payload": {
                        "content": QE_INPUT,
                        "mobile_indices": [0],
                        "use_environ": False,
                        "calc_mode": "continue",
                    },
                }
                self._send_json(200, payload)
                return

            if self.path == f"/api/zpe/compute/jobs/{self.job_id}/result":
                length = int(self.headers.get("Content-Length", "0"))
                body = self.rfile.read(length).decode("utf-8") if length else "{}"
                handler_cls.result_payload = json.loads(body)
                self._send_json(200, {"ok": True})
                return

            if self.path == f"/api/zpe/compute/jobs/{self.job_id}/failed":
                length = int(self.headers.get("Content-Length", "0"))
                body = self.rfile.read(length).decode("utf-8") if length else "{}"
                handler_cls.failed_payload = json.loads(body)
                self._send_json(200, {"ok": True})
                return

            self._send_json(404, {"detail": "not found"})

    return _LeaseHandler


def _start_server() -> tuple[HTTPServer, str, type[BaseHTTPRequestHandler]]:
    handler_cls = _make_handler()
    server = HTTPServer(("127.0.0.1", 0), handler_cls)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host_raw = server.server_address[0]
    port = server.server_address[1]
    host = host_raw.decode("utf-8") if isinstance(host_raw, bytes) else str(host_raw)
    return server, f"http://{host}:{port}", handler_cls


def test_http_worker_success(monkeypatch, tmp_path):
    server, base_url, handler_cls = _start_server()

    settings = ZPESettings(
        worker_mode="mock",
        work_dir=str(tmp_path),
    )
    monkeypatch.setattr(http_worker, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_worker, "get_zpe_settings", lambda: settings)

    try:
        lease = http_worker._lease_job(base_url, handler_cls.token, timeout=5)
        assert lease is not None
        artifacts = zpe_worker.compute_zpe_artifacts(
            lease["payload"], job_id=lease["job_id"]
        )
        http_worker._submit_result(
            base_url,
            handler_cls.token,
            lease["job_id"],
            "tenant-1",
            lease["lease_id"],
            artifacts.result,
            artifacts.summary_text,
            artifacts.freqs_csv,
            meta={"worker_hostname": "test"},
            execution_event={
                "event_id": f"execution:{lease['job_id']}:{lease['lease_id']}:completed",
                "tenant_id": "tenant-1",
                "workspace_id": "workspace-1",
                "job_id": lease["job_id"],
                "submission_id": "submission-1",
                "execution_id": "node-1:job-1:lease-1",
                "state": "completed",
                "occurred_at": "2026-01-01T00:00:00+00:00",
                "trace_id": "trace-req-1",
                "scheduler_ref": {
                    "slurm_job_id": "12345",
                    "partition": "short",
                    "qos": "normal",
                },
                "result_ref": {
                    "output_uri": f"zpe://jobs/{lease['job_id']}/result",
                    "metadata_uri": f"zpe://jobs/{lease['job_id']}/files/summary",
                },
            },
            timeout=5,
        )
    finally:
        server.shutdown()

    assert handler_cls.result_payload is not None
    assert handler_cls.result_payload["lease_id"] == handler_cls.lease_id
    # Keep management-node -> FastAPI payload compatibility explicit.
    validated = ComputeResultRequest(**handler_cls.result_payload)
    assert validated.tenant_id == "tenant-1"
    assert validated.lease_id == handler_cls.lease_id
    event = handler_cls.result_payload["execution_event"]
    assert event["event_id"] == f"execution:{handler_cls.job_id}:{handler_cls.lease_id}:completed"
    assert event["tenant_id"] == "tenant-1"
    assert event["workspace_id"] == "workspace-1"
    assert event["submission_id"] == "submission-1"
    assert event["state"] == "completed"
    assert event["trace_id"] == "trace-req-1"
    assert event["scheduler_ref"]["slurm_job_id"] == "12345"
    assert event["result_ref"]["output_uri"] == f"zpe://jobs/{handler_cls.job_id}/result"


def test_http_worker_failed(monkeypatch, tmp_path):
    server, base_url, handler_cls = _start_server()

    settings = ZPESettings(
        worker_mode="mock",
        work_dir=str(tmp_path),
    )
    monkeypatch.setattr(http_worker, "get_zpe_settings", lambda: settings)
    monkeypatch.setattr(zpe_worker, "get_zpe_settings", lambda: settings)

    try:
        lease = http_worker._lease_job(base_url, handler_cls.token, timeout=5)
        assert lease is not None
        http_worker._submit_failed(
            base_url,
            handler_cls.token,
            lease["job_id"],
            "tenant-1",
            lease["lease_id"],
            "TEST_ERROR",
            "boom",
            "trace",
            execution_event={
                "event_id": f"execution:{lease['job_id']}:{lease['lease_id']}:failed",
                "tenant_id": "tenant-1",
                "workspace_id": "workspace-1",
                "job_id": lease["job_id"],
                "submission_id": "submission-1",
                "execution_id": "node-1:job-1:lease-1",
                "state": "failed",
                "occurred_at": "2026-01-01T00:00:00+00:00",
                "trace_id": "trace-req-1",
                "error": {
                    "code": "COMPUTE_ERROR",
                    "message": "boom",
                    "retryable": True,
                },
            },
            timeout=5,
        )
    finally:
        server.shutdown()

    assert handler_cls.failed_payload is not None
    assert handler_cls.failed_payload["error_code"] == "TEST_ERROR"
    validated = ComputeFailedRequest(**handler_cls.failed_payload)
    assert validated.tenant_id == "tenant-1"
    assert validated.lease_id == handler_cls.lease_id
    event = handler_cls.failed_payload["execution_event"]
    assert event["event_id"] == f"execution:{handler_cls.job_id}:{handler_cls.lease_id}:failed"
    assert event["state"] == "failed"
    assert event["trace_id"] == "trace-req-1"
    assert event["error"]["code"] == "COMPUTE_ERROR"
    assert event["error"]["retryable"] is True
