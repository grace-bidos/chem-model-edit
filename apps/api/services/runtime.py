from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any, cast
from urllib import error as urlerror
from urllib import parse as urlparse
from urllib import request as urlrequest

from app.schemas.runtime import (
    ExecutionEvent,
    RuntimeJobStatusResponse,
    SubmitJobAccepted,
    SubmitJobCommand,
)
from app.schemas.zpe import ZPEJobStatus

from .runtime_settings import get_runtime_settings

class RuntimeConflictError(ValueError):
    pass


class RuntimeNotFoundError(KeyError):
    pass


class RuntimeConfigurationError(RuntimeError):
    pass


class RuntimeDownstreamError(RuntimeError):
    pass


@dataclass(frozen=True)
class SubmitResult:
    response: SubmitJobAccepted
    idempotent: bool


@dataclass(frozen=True)
class EventAck:
    idempotent: bool


def _require_url(value: str | None, *, key: str) -> str:
    if value and value.strip():
        return value.strip()
    raise RuntimeConfigurationError(f"missing runtime setting: {key}")


def _build_url(template: str, *, job_id: str) -> str:
    quoted = urlparse.quote(job_id, safe="")
    return template.format(job_id=quoted)


def _parse_json_body(payload: bytes | None) -> dict[str, Any]:
    if not payload:
        return {}
    text = payload.decode("utf-8", errors="replace").strip()
    if not text:
        return {}
    parsed_obj: object = json.loads(text)
    if not isinstance(parsed_obj, dict):
        return {}
    parsed_dict: dict[str, Any] = {}
    for key, value in parsed_obj.items():  # pyright: ignore[reportUnknownVariableType]
        if isinstance(key, str):
            parsed_dict[key] = value
    return parsed_dict


def _error_detail(payload: dict[str, Any]) -> str:
    if not payload:
        return "downstream_error"
    error = payload.get("error")
    if isinstance(error, dict):
        error_map = cast(dict[str, Any], error)
        message = error_map.get("message", None)
        if isinstance(message, str) and message.strip():
            return message.strip()
    detail = payload.get("detail")
    if isinstance(detail, str) and detail.strip():
        return detail.strip()
    return "downstream_error"


class RuntimeGateway:
    def __init__(self) -> None:
        self._settings = get_runtime_settings()
        timeout = self._settings.request_timeout_seconds
        self._timeout_seconds = timeout if timeout > 0 else 5

    def _headers(self, *, tenant_id: str) -> dict[str, str]:
        headers = {
            "Content-Type": "application/json",
            "x-tenant-id": tenant_id,
        }
        token = self._settings.service_auth_bearer_token
        if token:
            headers["Authorization"] = f"Bearer {token}"
        return headers

    def _request_json(
        self,
        *,
        method: str,
        url: str,
        tenant_id: str,
        body: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        encoded = None
        if body is not None:
            encoded = json.dumps(body, ensure_ascii=True, separators=(",", ":")).encode("utf-8")
        req = urlrequest.Request(
            url,
            method=method,
            data=encoded,
            headers=self._headers(tenant_id=tenant_id),
        )
        try:
            with urlrequest.urlopen(req, timeout=self._timeout_seconds) as resp:
                return _parse_json_body(resp.read())
        except urlerror.HTTPError as exc:
            payload = _parse_json_body(exc.read())
            detail = _error_detail(payload)
            if exc.code == 404:
                raise RuntimeNotFoundError(detail) from exc
            if exc.code in {409, 422}:
                raise RuntimeConflictError(detail) from exc
            raise RuntimeDownstreamError(
                f"runtime downstream HTTP {exc.code}: {detail}"
            ) from exc
        except urlerror.URLError as exc:
            raise RuntimeDownstreamError("runtime downstream request failed") from exc

    def submit(self, command: SubmitJobCommand, *, trace_id: str) -> SubmitResult:
        url = _require_url(
            self._settings.command_submit_url,
            key="RUNTIME_COMMAND_SUBMIT_URL",
        )
        payload = command.model_dump()
        payload["trace_id"] = trace_id
        response = self._request_json(
            method="POST",
            url=url,
            tenant_id=command.tenant_id,
            body=payload,
        )
        accepted = SubmitJobAccepted.model_validate(response)
        return SubmitResult(
            response=accepted,
            idempotent=bool(response.get("idempotent", False)),
        )

    def apply_event(self, event: ExecutionEvent) -> EventAck:
        template = _require_url(
            self._settings.command_event_url_template,
            key="RUNTIME_COMMAND_EVENT_URL_TEMPLATE",
        )
        url = _build_url(template, job_id=event.job_id)
        response = self._request_json(
            method="POST",
            url=url,
            tenant_id=event.tenant_id,
            body=event.model_dump(),
        )
        return EventAck(idempotent=bool(response.get("idempotent", False)))

    def get_status(self, job_id: str, *, tenant_id: str) -> RuntimeJobStatusResponse:
        template = _require_url(
            self._settings.read_status_url_template,
            key="RUNTIME_READ_STATUS_URL_TEMPLATE",
        )
        url = _build_url(template, job_id=job_id)
        payload = self._request_json(method="GET", url=url, tenant_id=tenant_id)
        return RuntimeJobStatusResponse.model_validate(payload)

    def get_detail(self, job_id: str, *, tenant_id: str) -> dict[str, Any]:
        template = _require_url(
            self._settings.read_detail_url_template,
            key="RUNTIME_READ_DETAIL_URL_TEMPLATE",
        )
        url = _build_url(template, job_id=job_id)
        return self._request_json(method="GET", url=url, tenant_id=tenant_id)

    def get_projection_status(self, job_id: str, *, tenant_id: str) -> ZPEJobStatus:
        template = _require_url(
            self._settings.read_projection_url_template,
            key="RUNTIME_READ_PROJECTION_URL_TEMPLATE",
        )
        url = _build_url(template, job_id=job_id)
        payload = self._request_json(method="GET", url=url, tenant_id=tenant_id)
        return ZPEJobStatus.model_validate(payload)


_runtime_gateway: RuntimeGateway | None = None


def get_runtime_store() -> RuntimeGateway:
    global _runtime_gateway
    if _runtime_gateway is None:
        _runtime_gateway = RuntimeGateway()
    return _runtime_gateway
