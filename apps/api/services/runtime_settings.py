from __future__ import annotations

from functools import lru_cache
from typing import Literal

from pydantic_settings import BaseSettings, SettingsConfigDict

from services.zpe.settings import resolve_env_file


class RuntimeSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="RUNTIME_",
        extra="ignore",
    )

    execution_owner: Literal["aiida-slurm", "mock"] = "aiida-slurm"
    request_timeout_seconds: int = 5
    service_auth_bearer_token: str | None = None

    command_submit_url: str | None = None
    command_event_url_template: str | None = None
    read_status_url_template: str | None = None
    read_detail_url_template: str | None = None
    read_projection_url_template: str | None = None


@lru_cache
def get_runtime_settings() -> RuntimeSettings:
    return RuntimeSettings(_env_file=resolve_env_file())  # pyright: ignore[reportCallIssue]
