from __future__ import annotations

from fastapi import APIRouter, HTTPException, Request

from app.deps import require_user_session
from app.schemas.auth import (
    AuthLoginRequest,
    AuthLogoutResponse,
    AuthMe,
    AuthRegisterRequest,
    AuthSession,
    AuthUser as AuthUserSchema,
)
from services import auth as auth_service
from services.auth.store import AuthUser

router = APIRouter(prefix="/api/auth", tags=["auth"])


def _to_auth_user(user: AuthUser) -> AuthUserSchema:
    """内部ユーザー型をAPIレスポンス型へ変換する．

    Parameters:
        user: 認証ストアのユーザー情報．

    Returns:
        APIレスポンス用のユーザー情報．
    """
    return AuthUserSchema(id=user.user_id, email=user.email, created_at=user.created_at)


@router.post("/register", response_model=AuthSession)
async def register(request: AuthRegisterRequest) -> AuthSession:
    store = auth_service.get_auth_store()
    try:
        user = store.create_user(request.email, request.password)
    except ValueError as exc:
        detail = str(exc)
        status = 409 if "already" in detail else 400
        raise HTTPException(status_code=status, detail=detail) from exc
    session = store.create_session(user.user_id)
    return AuthSession(
        token=session.token,
        expires_at=session.expires_at,
        user=_to_auth_user(user),
    )


@router.post("/login", response_model=AuthSession)
async def login(request: AuthLoginRequest) -> AuthSession:
    store = auth_service.get_auth_store()
    user = store.authenticate(request.email, request.password)
    if not user:
        raise HTTPException(status_code=401, detail="invalid credentials")
    session = store.create_session(user.user_id)
    return AuthSession(
        token=session.token,
        expires_at=session.expires_at,
        user=_to_auth_user(user),
    )


@router.post("/logout", response_model=AuthLogoutResponse)
async def logout(raw: Request) -> AuthLogoutResponse:
    store = auth_service.get_auth_store()
    token = raw.headers.get("authorization")
    if token and token.lower().startswith("bearer "):
        store.delete_session(token.split(" ", 1)[1].strip())
    return AuthLogoutResponse()


@router.get("/me", response_model=AuthMe)
async def me(raw: Request) -> AuthMe:
    user, session = require_user_session(raw)
    return AuthMe(user=_to_auth_user(user), expires_at=session.expires_at)
