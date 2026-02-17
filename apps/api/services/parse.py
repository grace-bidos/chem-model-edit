from __future__ import annotations

from importlib import import_module
from io import StringIO
import re
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, SupportsFloat, Tuple, cast

from ase import Atoms as ASEAtoms
from pymatgen.io.pwscf import PWInput

from app.schemas.common import Atom, Lattice, QeParameters, Structure, Vector3

ase_read = cast(Callable[..., Any], getattr(import_module("ase.io"), "read"))


def _as_float3(value: object) -> tuple[float, float, float]:
    vec = cast(tuple[object, object, object], value)
    return (
        float(cast(SupportsFloat, vec[0])),
        float(cast(SupportsFloat, vec[1])),
        float(cast(SupportsFloat, vec[2])),
    )


def _ase_symbols_positions(atoms_obj: ASEAtoms) -> list[tuple[str, tuple[float, float, float]]]:
    symbols = cast(list[str], cast(Any, atoms_obj).get_chemical_symbols())
    positions = cast(list[object], cast(Any, atoms_obj).get_positions())
    return [
        (symbol, _as_float3(position))
        for symbol, position in zip(symbols, positions, strict=True)
    ]


def _ase_lattice_vectors(atoms_obj: ASEAtoms) -> list[list[float]]:
    cell = cast(Any, atoms_obj).get_cell()
    raw_vectors = cast(list[list[object]], cell.tolist())
    return [
        [float(cast(SupportsFloat, row[0])), float(cast(SupportsFloat, row[1])), float(cast(SupportsFloat, row[2]))]
        for row in raw_vectors
    ]


def _ase_atoms_from_content(content: str) -> ASEAtoms:
    atoms_result: object = ase_read(StringIO(content), format="espresso-in")
    if isinstance(atoms_result, list):
        if not atoms_result:
            raise ValueError("ASEが構造を返しませんでした。")
        return cast(ASEAtoms, atoms_result[0])
    return cast(ASEAtoms, atoms_result)


def _pymatgen_atoms_from_content(content: str) -> ASEAtoms:
    pw_input = PWInput.from_str(content)
    structure = cast(Any, pw_input).structure
    sites = cast(Sequence[Any], structure.sites)
    symbols = [
        getattr(site.specie, "symbol", str(site.specie))
        for site in sites
    ]
    positions = [
        (
            float(cast(SupportsFloat, site.coords[0])),
            float(cast(SupportsFloat, site.coords[1])),
            float(cast(SupportsFloat, site.coords[2])),
        )
        for site in sites
    ]
    cell = cast(list[list[float]], structure.lattice.matrix.tolist())
    return ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=True,
    )


def structure_from_ase(atoms_obj: ASEAtoms) -> Structure:
    parsed: List[Atom] = []
    for symbol, (x, y, z) in _ase_symbols_positions(atoms_obj):
        parsed.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
    lattice = _lattice_from_vectors(_ase_lattice_vectors(atoms_obj))
    return Structure(atoms=parsed, lattice=lattice)


def _from_ase(content: str) -> Structure:
    atoms_obj = _ase_atoms_from_content(content)
    return structure_from_ase(atoms_obj)


def _from_pymatgen(content: str) -> Structure:
    atoms_obj = _pymatgen_atoms_from_content(content)
    return structure_from_ase(atoms_obj)


_POSITION_ANY = re.compile(r"^\s*ATOMIC_POSITIONS\b", re.IGNORECASE | re.MULTILINE)


def _lattice_from_vectors(
    vectors: Optional[Sequence[Sequence[float]]],
) -> Optional[Lattice]:
    if not vectors or len(vectors) < 3:
        return None
    flat = [abs(value) for row in vectors[:3] for value in row[:3]]
    if all(value < 1.0e-8 for value in flat):
        return None
    a, b, c = vectors[0], vectors[1], vectors[2]
    return Lattice(
        a=Vector3(x=float(a[0]), y=float(a[1]), z=float(a[2])),
        b=Vector3(x=float(b[0]), y=float(b[1]), z=float(b[2])),
        c=Vector3(x=float(c[0]), y=float(c[1]), z=float(c[2])),
    )


def parse_qe_in(content: str) -> Structure:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        return _from_ase(content)
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _from_pymatgen(content)
        except Exception as pmg_error:
            raise ValueError(
                f"QE .in のパースに失敗しました: ASE={ase_error}, pymatgen={pmg_error}"
            ) from pmg_error


def parse_qe_atoms(content: str) -> Tuple[ASEAtoms, str]:
    if not _POSITION_ANY.search(content):
        raise ValueError("構造データではありません。")
    try:
        return _ase_atoms_from_content(content), "ase"
    except Exception as ase_error:  # pragma: no cover - fallback path varies
        try:
            return _pymatgen_atoms_from_content(content), "pymatgen"
        except Exception as pmg_error:
            raise ValueError(
                f"QE .in のパースに失敗しました: ASE={ase_error}, pymatgen={pmg_error}"
            ) from pmg_error


def _safe_mapping(value: object) -> Dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, Mapping):
        mapping = cast(Mapping[object, Any], value)
        return {str(key): item for key, item in mapping.items()}
    try:
        mapped = cast(Mapping[object, Any], dict(cast(Any, value)))
    except Exception:
        return {}
    return {str(key): item for key, item in mapped.items()}


def _safe_dict_like(value: object) -> Optional[Dict[str, Any]]:
    if value is None:
        return None
    for attr in ("as_dict", "to_dict"):
        fn = getattr(value, attr, None)
        if callable(fn):
            try:
                result = fn()
            except Exception:
                return None
            if isinstance(result, Mapping):
                mapping = cast(Mapping[object, Any], result)
                return {str(key): item for key, item in mapping.items()}
            return None
    return None


def extract_qe_params(content: str) -> Optional[QeParameters]:
    try:
        pw_input = PWInput.from_str(content)
    except Exception:
        return None

    control = _safe_mapping(getattr(pw_input, "control", None))
    system = _safe_mapping(getattr(pw_input, "system", None))
    electrons = _safe_mapping(getattr(pw_input, "electrons", None))
    ions = _safe_mapping(getattr(pw_input, "ions", None))
    cell = _safe_mapping(getattr(pw_input, "cell", None))

    pseudos_raw = getattr(pw_input, "pseudopotentials", None)
    if pseudos_raw is None:
        pseudos_raw = getattr(pw_input, "pseudo", None)
    pseudos = {
        str(key): str(value) for key, value in _safe_mapping(pseudos_raw).items()
    }

    kpoints = _safe_dict_like(getattr(pw_input, "kpoints", None))

    return QeParameters(
        control=control,
        system=system,
        electrons=electrons,
        ions=ions,
        cell=cell,
        pseudopotentials=pseudos,
        kpoints=kpoints,
    )
