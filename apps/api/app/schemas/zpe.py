from __future__ import annotations

from typing import Any, Dict, List, Literal, Optional, Tuple

from pydantic import Field

from .base import ApiModel
from .common import Pagination, Structure


JobState = Literal["queued", "started", "finished", "failed"]


class ZPEParseRequest(ApiModel):
    content: str
    structure_id: Optional[str] = None


class ZPEParseResponse(ApiModel):
    structure: Structure
    fixed_indices: List[int]
    atomic_species: Dict[str, str] = Field(default_factory=dict)
    kpoints: Optional[Tuple[int, int, int]] = None


class ZPEJobRequest(ApiModel):
    calc_type: Literal["qe.zpe.v1", "qe.relax.v1"] = "qe.zpe.v1"
    content: str
    mobile_indices: List[int]
    use_environ: bool = False
    input_dir: Optional[str] = None
    calc_mode: Literal["new", "continue"] = "continue"
    structure_id: Optional[str] = None


class ZPEJobResponse(ApiModel):
    id: str


class ZPEJobStatus(ApiModel):
    status: JobState
    detail: Optional[str] = None
    updated_at: Optional[str] = None


class ZPEResult(ApiModel):
    calc_type: Literal["qe.zpe.v1", "qe.relax.v1"] = "qe.zpe.v1"
    freqs_cm: List[float]
    zpe_ev: float
    s_vib_jmol_k: float
    mobile_indices: List[int]
    fixed_indices: List[int]
    kpoints: Tuple[int, int, int]
    delta: float
    low_cut_cm: float
    temperature: float
    use_environ: bool
    qe_input: str
    pseudo_dir: str
    calc_start_time: str
    calc_end_time: str
    elapsed_seconds: float
    cache_checked: int
    cache_deleted: int
    ecutwfc: Optional[float] = None
    ecutrho: Optional[float] = None


class ZPEJobResultResponse(ApiModel):
    result: ZPEResult


class EnrollTokenRequest(ApiModel):
    ttl_seconds: Optional[int] = None
    label: Optional[str] = None


class EnrollTokenResponse(ApiModel):
    token: str
    expires_at: str
    ttl_seconds: int
    label: Optional[str] = None


class ComputeRegisterRequest(ApiModel):
    token: str
    name: Optional[str] = None
    queue_name: Optional[str] = None
    meta: Dict[str, Any] = Field(default_factory=dict)


class ComputeRegisterResponse(ApiModel):
    id: str
    registered_at: str
    name: Optional[str] = None
    worker_token: str
    token_expires_at: str
    token_ttl_seconds: int


class ComputeRevokeResponse(ApiModel):
    revoked_count: int


class QueueTarget(ApiModel):
    id: str
    queue_name: str
    server_id: str
    registered_at: str
    name: Optional[str] = None


class QueueTargetListResponse(ApiModel):
    targets: List[QueueTarget]
    active_target_id: Optional[str] = None
    pagination: Pagination


class QueueTargetSelectResponse(ApiModel):
    active_target_id: str


class ComputeLeaseResponse(ApiModel):
    job_id: str
    payload: Dict[str, Any]
    lease_id: str
    lease_ttl_seconds: int
    meta: Dict[str, Any] = Field(default_factory=dict)


class ComputeResultRequest(ApiModel):
    lease_id: str
    result: Dict[str, Any]
    summary_text: str
    freqs_csv: str
    meta: Dict[str, Any] = Field(default_factory=dict)


class ComputeResultResponse(ApiModel):
    ok: bool = True
    idempotent: bool = False


class ComputeFailedRequest(ApiModel):
    lease_id: str
    error_code: str
    error_message: str
    traceback: Optional[str] = None


class ComputeFailedResponse(ApiModel):
    ok: bool = True
    requeued: bool
    retry_count: int


class OpsFlagsRequest(ApiModel):
    submission_enabled: Optional[bool] = None
    dequeue_enabled: Optional[bool] = None
    submission_route: Optional[Literal["redis-worker", "next-gen"]] = None
    result_read_source: Optional[Literal["redis", "projection"]] = None
    legacy_worker_endpoints_enabled: Optional[bool] = None


class OpsFlagsResponse(ApiModel):
    submission_enabled: bool
    dequeue_enabled: bool
    submission_route: Literal["redis-worker", "next-gen"]
    result_read_source: Literal["redis", "projection"]
    legacy_worker_endpoints_enabled: bool
