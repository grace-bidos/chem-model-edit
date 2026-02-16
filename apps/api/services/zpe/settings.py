from __future__ import annotations

from functools import lru_cache
import os
from pathlib import Path
from typing import Literal, Optional

from pydantic_settings import BaseSettings, SettingsConfigDict


def resolve_env_file() -> str:
    override = os.getenv("ZPE_ENV_FILE")
    if override:
        return override
    cwd = Path.cwd()
    for base in (cwd, *cwd.parents):
        candidate = base / ".env"
        if candidate.is_file():
            return str(candidate)
    return ".env"


class ZPESettings(BaseSettings):
    model_config = SettingsConfigDict(
        env_prefix="ZPE_",
        env_file=resolve_env_file(),
        extra="ignore",
    )

    redis_url: str = "redis://localhost:6379/0"
    work_dir: str = "zpe_jobs"
    queue_name: str = "zpe"
    compute_mode: str = "remote-queue"
    result_store: str = "redis"
    worker_mode: Literal["qe", "mock"] = "qe"
    admin_token: Optional[str] = None
    enroll_token_ttl_seconds: int = 3600
    worker_token_ttl_seconds: int = 604800
    submission_enabled: bool = True
    dequeue_enabled: bool = True
    cutover_submission_route: Literal["redis-worker", "next-gen"] = "redis-worker"
    cutover_result_read_source: Literal["redis", "projection"] = "redis"
    legacy_worker_endpoints_enabled: bool = True
    ops_flag_ttl_seconds: int = 3600

    pw_command: str = "pw.x"
    pw_path: Optional[str] = None
    use_mpi: bool = True
    mpi_cmd: str = "mpirun"
    np_core: int = 12

    pseudo_dir: Optional[str] = None
    allow_input_pseudo_dir: bool = False

    delta: float = 0.01
    low_cut_cm: float = 50.0
    temperature: float = 298.15
    vib_name: str = "vib_ads"
    outdir_name: str = "calc_scratch"
    calc_dir_name: str = "vib_ads_vac"
    job_timeout_seconds: int = 86400
    result_ttl_seconds: int = 604800
    lease_ttl_seconds: int = 600
    retry_max: int = 3
    retry_base_delay_seconds: int = 10
    retry_max_delay_seconds: int = 300

    environ_path: Optional[str] = None
    control_api_url: str = "http://localhost:8000"
    worker_token: Optional[str] = None
    worker_poll_interval_seconds: int = 10
    worker_poll_max_interval_seconds: int = 60
    worker_request_timeout_seconds: int = 60
    convex_relay_url: Optional[str] = None
    convex_relay_token: Optional[str] = None
    convex_relay_timeout_seconds: int = 5
    slurm_policy_path: Optional[str] = None
    slurm_adapter: Literal["stub-policy", "passthrough", "real-policy"] = "stub-policy"
    slurm_adapter_rollback_guard: Literal[
        "allow",
        "force-stub-policy",
        "force-passthrough",
    ] = "allow"
    slurm_real_adapter_probe_timeout_seconds: int = 5


@lru_cache
def get_zpe_settings() -> ZPESettings:
    """ZPE設定を環境変数から読み込み，キャッシュして返す．

    Returns:
        解決済みの設定オブジェクト．
    """
    return ZPESettings()
