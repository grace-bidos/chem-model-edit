from __future__ import annotations

import secrets
from typing import Tuple

from fastapi import HTTPException, Request

from services import auth as auth_service
from services.authn import UserIdentity, get_authn_settings, verify_clerk_token
from services.auth.store import AuthSession, AuthUser
from services.zpe.job_owner import get_job_owner_store
from services.zpe.settings import get_zpe_settings
from services.zpe.worker_auth import get_worker_token_store


def _extract_bearer_token(request: Request) -> str | None:
    auth = request.headers.get("authorization")
    if auth and auth.lower().startswith("bearer "):
        token = auth.split(" ", 1)[1].strip()
        return token or None
    return None


def require_user_session(request: Request) -> Tuple[AuthUser, AuthSession]:
    store = auth_service.get_auth_store()
    token = _extract_bearer_token(request)
    if not token:
        raise HTTPException(status_code=401, detail="unauthorized")
    session = store.get_user_by_session(token, refresh=True)
    if not session:
        raise HTTPException(status_code=401, detail="unauthorized")
    user = store.get_user_by_id(session.user_id)
    if not user:
        raise HTTPException(status_code=401, detail="unauthorized")
    return user, session


def _require_dev_bypass_identity(request: Request) -> UserIdentity:
    settings = get_authn_settings()
    user_id = request.headers.get(settings.dev_bypass_user_id_header)
    if not user_id:
        raise HTTPException(status_code=401, detail="unauthorized")
    email = request.headers.get(settings.dev_bypass_email_header)
    allowed = settings.allowed_email_set()
    if allowed and (not email or email.lower() not in allowed):
        raise HTTPException(status_code=403, detail="forbidden")
    return UserIdentity(user_id=user_id, email=email)


def require_user_identity(request: Request) -> UserIdentity:
    mode = get_authn_settings().mode
    if mode == "dev-bypass":
        return _require_dev_bypass_identity(request)
    if mode == "clerk":
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
    user, _session = require_user_session(request)
    return UserIdentity(user_id=user.user_id, email=user.email)


def require_job_owner(request: Request, job_id: str) -> UserIdentity:
    user = require_user_identity(request)
    owner_store = get_job_owner_store()
    owner_id = owner_store.get_owner(job_id)
    if not owner_id:
        raise HTTPException(status_code=404, detail="job not found")
    if owner_id != user.user_id:
        raise HTTPException(status_code=403, detail="forbidden")
    return user


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


def require_worker(request: Request) -> str:
    token = _extract_bearer_token(request)
    if not token:
        raise HTTPException(status_code=401, detail="missing worker token")
    store = get_worker_token_store()
    try:
        return store.validate(token)
    except PermissionError as exc:
        detail = str(exc) or "invalid token"
        status = 403 if "revoked" in detail else 401
        raise HTTPException(status_code=status, detail=detail) from exc
