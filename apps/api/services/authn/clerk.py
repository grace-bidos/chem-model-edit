from __future__ import annotations

from functools import lru_cache
from typing import Any

import jwt
from jwt import InvalidTokenError, PyJWKClient

from .settings import get_authn_settings
from .types import UserIdentity


def _jwks_url(issuer: str) -> str:
    return f"{issuer.rstrip('/')}/.well-known/jwks.json"


@lru_cache(maxsize=8)
def _get_jwks_client(issuer: str) -> PyJWKClient:
    return PyJWKClient(_jwks_url(issuer))


def _decode_clerk_jwt(
    token: str,
    *,
    issuer: str,
    audience: str | None,
) -> dict[str, Any]:
    client = _get_jwks_client(issuer)
    signing_key = client.get_signing_key_from_jwt(token)
    decode_kwargs: dict[str, Any] = {
        "key": signing_key.key,
        "algorithms": ["RS256", "RS384", "RS512", "EdDSA"],
        "issuer": issuer,
        "options": {"verify_aud": bool(audience)},
    }
    if audience:
        decode_kwargs["audience"] = audience
    return jwt.decode(token, **decode_kwargs)


def verify_clerk_token(token: str) -> UserIdentity:
    settings = get_authn_settings()
    issuer = settings.clerk_issuer
    if not issuer:
        raise PermissionError("clerk issuer not configured")

    try:
        claims = _decode_clerk_jwt(
            token,
            issuer=issuer,
            audience=settings.clerk_audience,
        )
    except InvalidTokenError as exc:
        raise PermissionError("unauthorized") from exc

    user_id = claims.get("sub")
    if not isinstance(user_id, str) or not user_id:
        raise PermissionError("unauthorized")

    email_value = claims.get("email")
    email = email_value.lower() if isinstance(email_value, str) else None

    allowed = settings.allowed_email_set()
    if allowed and (not email or email not in allowed):
        raise PermissionError("allowlist denied")

    return UserIdentity(user_id=user_id, email=email)
