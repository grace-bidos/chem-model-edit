from __future__ import annotations

from functools import lru_cache
from typing import Literal

from pydantic_settings import BaseSettings, SettingsConfigDict


class AuthnSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="AUTH_",
        env_file=".env",
        extra="ignore",
    )

    mode: Literal["clerk", "dev-bypass"] = "clerk"
    clerk_issuer: str | None = None
    clerk_audience: str | None = None
    allowed_emails: str | None = None
    dev_bypass_user_id_header: str = "X-Dev-User-Id"
    dev_bypass_email_header: str = "X-Dev-User-Email"

    def allowed_email_set(self) -> set[str]:
        if not self.allowed_emails:
            return set()
        return {
            token.strip().lower()
            for token in self.allowed_emails.split(",")
            if token.strip()
        }


@lru_cache(maxsize=1)
def get_authn_settings() -> AuthnSettings:
    return AuthnSettings()
