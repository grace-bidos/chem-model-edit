from __future__ import annotations

import pytest

from services.auth.passwords import _parse_hash, hash_password, verify_password


def test_hash_format_parses_and_verifies() -> None:
    stored = hash_password("correct-horse-battery-staple", iterations=1_500)

    iterations, salt, digest = _parse_hash(stored)

    assert iterations == 1_500
    assert len(salt) == 16
    assert len(digest) == 32
    assert verify_password("correct-horse-battery-staple", stored) is True
    assert verify_password("wrong-password", stored) is False


def test_verify_requires_matching_pepper() -> None:
    stored = hash_password("swordfish-pass", iterations=1_200, pepper="pepper-A")

    assert verify_password("swordfish-pass", stored, pepper="pepper-A") is True
    assert verify_password("swordfish-pass", stored, pepper="pepper-B") is False
    assert verify_password("swordfish-pass", stored) is False


@pytest.mark.parametrize(
    ("stored", "message"),
    [
        ("not-a-split-hash", "invalid password hash format"),
        (
            "sha256$1200$c2FsdA==$ZGlnZXN0",
            "unsupported password hash algorithm",
        ),
    ],
)
def test_parse_hash_invalid_format_handling(stored: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        _parse_hash(stored)


def test_verify_invalid_hash_format_raises() -> None:
    with pytest.raises(ValueError, match="invalid password hash format"):
        verify_password("pw", "broken-hash")
