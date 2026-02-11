from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class UserIdentity:
    user_id: str
    email: str | None = None
