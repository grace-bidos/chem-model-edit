from __future__ import annotations

from types import SimpleNamespace

import fakeredis
from hypothesis import assume, given
from hypothesis import strategies as st
import pytest

import services.zpe.submit_idempotency as submit_idempotency
from services.zpe.submit_idempotency import SubmitIdempotencyRecord, SubmitIdempotencyStore

_ID_RE = r"[a-z][a-z0-9_-]{0,12}"
_FP_RE = r"[a-f0-9]{64}"


@pytest.fixture(autouse=True)
def _patch_zpe_settings(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        submit_idempotency,
        "get_zpe_settings",
        lambda: SimpleNamespace(result_ttl_seconds=600),
    )


@given(
    tenant_id=st.from_regex(_ID_RE, fullmatch=True),
    user_id=st.from_regex(_ID_RE, fullmatch=True),
    request_id=st.from_regex(_ID_RE, fullmatch=True),
    request_fingerprint=st.from_regex(_FP_RE, fullmatch=True),
    job_id=st.from_regex(_ID_RE, fullmatch=True),
    requested_queue_name=st.from_regex(_ID_RE, fullmatch=True),
    resolved_queue_name=st.from_regex(_ID_RE, fullmatch=True),
)
def test_claim_finalize_round_trip_returns_ready_record(
    tenant_id: str,
    user_id: str,
    request_id: str,
    request_fingerprint: str,
    job_id: str,
    requested_queue_name: str,
    resolved_queue_name: str,
) -> None:
    store = SubmitIdempotencyStore(redis=fakeredis.FakeRedis())

    claim = store.claim_or_get(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        request_fingerprint=request_fingerprint,
    )
    assert claim.state == "claimed"
    assert claim.claim_token is not None

    record = SubmitIdempotencyRecord(
        job_id=job_id,
        request_fingerprint=request_fingerprint,
        requested_queue_name=requested_queue_name,
        resolved_queue_name=resolved_queue_name,
    )
    store.finalize_claim(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        claim_token=claim.claim_token,
        record=record,
    )

    claimed_again = store.claim_or_get(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        request_fingerprint=request_fingerprint,
    )
    assert claimed_again.state == "ready"
    assert claimed_again.record == record


@given(
    tenant_id=st.from_regex(_ID_RE, fullmatch=True),
    user_id=st.from_regex(_ID_RE, fullmatch=True),
    request_id=st.from_regex(_ID_RE, fullmatch=True),
    first_fingerprint=st.from_regex(_FP_RE, fullmatch=True),
    second_fingerprint=st.from_regex(_FP_RE, fullmatch=True),
)
def test_claim_or_get_rejects_fingerprint_conflicts_for_same_request_id(
    tenant_id: str,
    user_id: str,
    request_id: str,
    first_fingerprint: str,
    second_fingerprint: str,
) -> None:
    assume(first_fingerprint != second_fingerprint)
    store = SubmitIdempotencyStore(redis=fakeredis.FakeRedis())

    store.claim_or_get(
        user_id=user_id,
        tenant_id=tenant_id,
        request_id=request_id,
        request_fingerprint=first_fingerprint,
    )

    with pytest.raises(ValueError, match="different submission payload"):
        store.claim_or_get(
            user_id=user_id,
            tenant_id=tenant_id,
            request_id=request_id,
            request_fingerprint=second_fingerprint,
        )
