from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic_settings import BaseSettings, SettingsConfigDict

from services.zpe.settings import resolve_env_file


class RuntimeSettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="RUNTIME_",
        extra="ignore",
    )

    sqlite_path: str = ".just-runtime/runtime_jobs.sqlite3"
    execution_owner: Literal["aiida-slurm", "mock"] = "aiida-slurm"


@lru_cache
def get_runtime_settings() -> RuntimeSettings:
    return RuntimeSettings(_env_file=resolve_env_file())  # pyright: ignore[reportCallIssue]


def resolve_runtime_db_path() -> Path:
    path = Path(get_runtime_settings().sqlite_path)
    if path.is_absolute():
        return path
    return Path.cwd() / path
