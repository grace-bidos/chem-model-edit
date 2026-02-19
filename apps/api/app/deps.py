from __future__ import annotations

import re
import secrets

from fastapi import HTTPException, Request

from services.authn import UserIdentity, get_authn_settings, verify_clerk_token
from services.zpe.settings import get_zpe_settings

_TENANT_HEADER = "x-tenant-id"
_TENANT_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{1,127}$")


def _extract_bearer_token(request: Request) -> str | None:
    auth = request.headers.get("authorization")
    if auth and auth.lower().startswith("bearer "):
        token = auth.split(" ", 1)[1].strip()
        return token or None
    return None


def _extract_tenant_id(request: Request) -> str | None:
    state_value = getattr(request.state, "tenant_id", None)
    if isinstance(state_value, str):
        value = state_value.strip()
        return value or None
    raw = request.headers.get(_TENANT_HEADER)
    if raw is None:
        return None
    value = raw.strip()
    return value or None


def require_tenant_id(request: Request) -> str:
    tenant_id = _extract_tenant_id(request)
    if tenant_id is None:
        raise HTTPException(status_code=400, detail="missing tenant_id")
    if not _TENANT_ID_PATTERN.fullmatch(tenant_id):
        raise HTTPException(status_code=400, detail="invalid tenant_id")
    request.state.tenant_id = tenant_id
    return tenant_id


def require_user_identity(request: Request) -> UserIdentity:
    require_tenant_id(request)
    _ = get_authn_settings()
    token = _extract_bearer_token(request)
    if not token:
        raise HTTPException(status_code=401, detail="unauthorized")
    try:
        return verify_clerk_token(token)
    except PermissionError as exc:
        detail = str(exc) or "unauthorized"
        status = 403 if "allowlist" in detail else 401
        raise HTTPException(
            status_code=status,
            detail="forbidden" if status == 403 else "unauthorized",
        ) from exc


def require_admin(request: Request) -> None:
    settings = get_zpe_settings()
    if not settings.admin_token:
        return
    token = _extract_bearer_token(request)
    if not token or not secrets.compare_digest(token, settings.admin_token):
        raise HTTPException(status_code=401, detail="unauthorized")


def has_admin_access(request: Request) -> bool:
    settings = get_zpe_settings()
    if not settings.admin_token:
        return False
    token = _extract_bearer_token(request)
    return bool(token and secrets.compare_digest(token, settings.admin_token))

