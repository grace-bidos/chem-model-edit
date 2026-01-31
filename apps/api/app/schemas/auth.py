from __future__ import annotations

from .base import ApiModel


class AuthRegisterRequest(ApiModel):
    email: str
    password: str


class AuthLoginRequest(ApiModel):
    email: str
    password: str


class AuthUser(ApiModel):
    id: str
    email: str
    created_at: str


class AuthSession(ApiModel):
    token: str
    expires_at: str
    user: AuthUser


class AuthMe(ApiModel):
    user: AuthUser
    expires_at: str


class AuthLogoutResponse(ApiModel):
    ok: bool = True
