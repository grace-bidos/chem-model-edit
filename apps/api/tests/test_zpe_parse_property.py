from __future__ import annotations

import re

from hypothesis import given
from hypothesis import strategies as st
import pytest

from services.zpe.parse import ensure_mobile_indices, parse_kpoints_automatic, parse_namelist


_KEY_RE = re.compile(r"^[a-z][a-z0-9_]{0,8}$")


@given(
    natoms=st.integers(min_value=1, max_value=32),
    mobile_raw=st.sets(st.integers(min_value=0, max_value=31), min_size=1, max_size=20),
    fixed_raw=st.sets(st.integers(min_value=0, max_value=31), max_size=20),
)
def test_ensure_mobile_indices_returns_sorted_unique_disjoint(
    natoms: int, mobile_raw: set[int], fixed_raw: set[int]
) -> None:
    mobile = [idx for idx in mobile_raw if idx < natoms]
    fixed = [idx for idx in fixed_raw if idx < natoms and idx not in mobile]
    if not mobile:
        mobile = [0]
        if 0 in fixed:
            fixed.remove(0)

    result = ensure_mobile_indices(mobile, natoms=natoms, fixed_indices=fixed)
    assert result == sorted(set(mobile))
    assert not set(result).intersection(fixed)


@given(
    natoms=st.integers(min_value=1, max_value=32),
    out_of_range=st.one_of(
        st.integers(max_value=-1),
        st.integers(min_value=33, max_value=200),
    ),
)
def test_ensure_mobile_indices_rejects_out_of_range(natoms: int, out_of_range: int) -> None:
    with pytest.raises(ValueError, match="範囲外"):
        ensure_mobile_indices([0, out_of_range], natoms=natoms, fixed_indices=[])


@given(
    entries=st.dictionaries(
        keys=st.from_regex(_KEY_RE, fullmatch=True),
        values=st.integers(min_value=-10_000, max_value=10_000),
        min_size=1,
        max_size=8,
    )
)
def test_parse_namelist_parses_integer_entries(entries: dict[str, int]) -> None:
    body = "\n".join(f"  {k} = {v}," for k, v in entries.items())
    content = f"&SYSTEM\n{body}\n/\n"
    parsed = parse_namelist(content, "system")
    assert parsed == {k.lower(): v for k, v in entries.items()}


@given(
    nkx=st.integers(min_value=-20, max_value=20),
    nky=st.integers(min_value=-20, max_value=20),
    nkz=st.integers(min_value=-20, max_value=20),
    sx=st.integers(min_value=-4, max_value=4),
    sy=st.integers(min_value=-4, max_value=4),
    sz=st.integers(min_value=-4, max_value=4),
)
def test_parse_kpoints_automatic_round_trip(
    nkx: int, nky: int, nkz: int, sx: int, sy: int, sz: int
) -> None:
    content = (
        "&CONTROL\n  calculation='scf'\n/\n"
        "K_POINTS (automatic)\n"
        f"{nkx} {nky} {nkz} {sx} {sy} {sz}\n"
    )
    parsed = parse_kpoints_automatic(content)
    assert parsed == ((nkx, nky, nkz), (sx, sy, sz))
