from __future__ import annotations

import re

from hypothesis import assume, given
from hypothesis import strategies as st
import pytest

from services.zpe.parse import ensure_mobile_indices, parse_kpoints_automatic, parse_namelist


_KEY_RE = re.compile(r"^[a-z][a-z0-9_]{0,8}$")
_SAFE_STRING_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,10}$")


def _render_namelist_value(value: object, uppercase_bool: bool) -> str:
    if isinstance(value, bool):
        if uppercase_bool:
            return ".TRUE." if value else ".FALSE."
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        # Use Fortran-style exponents to assert D-notation handling.
        return f"{value:.1f}d0"
    if isinstance(value, str):
        return f"'{value}'"
    raise TypeError(f"Unsupported value type: {type(value)!r}")


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
    mobile_raw=st.lists(st.integers(min_value=0, max_value=31), min_size=1, max_size=40),
    fixed_raw=st.lists(st.integers(min_value=0, max_value=31), max_size=40),
)
def test_ensure_mobile_indices_is_order_invariant_and_idempotent(
    natoms: int, mobile_raw: list[int], fixed_raw: list[int]
) -> None:
    mobile = [idx for idx in mobile_raw if idx < natoms]
    assume(bool(mobile))
    fixed = [idx for idx in fixed_raw if idx < natoms and idx not in mobile]
    canonical = ensure_mobile_indices(mobile, natoms=natoms, fixed_indices=fixed)
    permuted = ensure_mobile_indices(list(reversed(mobile)) + mobile[:2], natoms=natoms, fixed_indices=fixed)
    assert permuted == canonical
    assert ensure_mobile_indices(canonical, natoms=natoms, fixed_indices=fixed) == canonical


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
    natoms=st.integers(min_value=1, max_value=32),
    overlap=st.integers(min_value=0, max_value=31),
)
def test_ensure_mobile_indices_rejects_overlap_with_fixed(natoms: int, overlap: int) -> None:
    assume(overlap < natoms)
    with pytest.raises(ValueError, match="固定原子"):
        ensure_mobile_indices([overlap], natoms=natoms, fixed_indices=[overlap])


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
    entries=st.lists(
        st.tuples(
            st.from_regex(_KEY_RE, fullmatch=True),
            st.one_of(
                st.booleans(),
                st.integers(min_value=-10_000, max_value=10_000),
                st.integers(min_value=-999, max_value=999)
                .filter(lambda n: n % 10 != 0)
                .map(lambda n: n / 10.0),
                st.from_regex(_SAFE_STRING_RE, fullmatch=True),
            ),
        ),
        min_size=1,
        max_size=12,
    )
)
def test_parse_namelist_supports_scalar_types_and_last_assignment_wins(
    entries: list[tuple[str, object]]
) -> None:
    lines: list[str] = []
    expected: dict[str, object] = {}
    for i, (key, raw_value) in enumerate(entries):
        value = _render_namelist_value(raw_value, uppercase_bool=(i % 2 == 1))
        rendered_key = key.upper() if i % 2 == 1 else key
        lines.append(f"  {rendered_key} = {value}, ! inline comment")
        expected[key] = raw_value
    body = "\n".join(lines)
    content = f"&SYSTEM\n  ignored_chunk_without_equals,\n{body}\n/\n"
    parsed = parse_namelist(content, "SyStEm")
    assert parsed == expected


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


@given(
    nkx=st.integers(min_value=-20, max_value=20),
    nky=st.integers(min_value=-20, max_value=20),
    nkz=st.integers(min_value=-20, max_value=20),
    sx=st.integers(min_value=-4, max_value=4),
    sy=st.integers(min_value=-4, max_value=4),
    sz=st.integers(min_value=-4, max_value=4),
    header=st.sampled_from(["K_POINTS (automatic)", "k_points(Automatic)", "K_POINTS( AUTOMATIC )"]),
    spacer=st.sampled_from([" ", "  ", "\t"]),
)
def test_parse_kpoints_automatic_is_case_and_spacing_tolerant(
    nkx: int, nky: int, nkz: int, sx: int, sy: int, sz: int, header: str, spacer: str
) -> None:
    content = (
        "&CONTROL\n  calculation='scf'\n/\n"
        f"{header}\n"
        f"{nkx}{spacer}{nky}{spacer}{nkz}{spacer}{sx}{spacer}{sy}{spacer}{sz}\n"
    )
    assert parse_kpoints_automatic(content) == ((nkx, nky, nkz), (sx, sy, sz))


@given(
    mode=st.sampled_from(["gamma", "tpiba", "crystal", "tpiba_b"]),
    nkx=st.integers(min_value=1, max_value=20),
    nky=st.integers(min_value=1, max_value=20),
    nkz=st.integers(min_value=1, max_value=20),
)
def test_parse_kpoints_automatic_returns_none_for_non_automatic_modes(
    mode: str, nkx: int, nky: int, nkz: int
) -> None:
    content = f"K_POINTS ({mode})\n{nkx} {nky} {nkz}\n"
    assert parse_kpoints_automatic(content) is None
