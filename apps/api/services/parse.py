from __future__ import annotations

from io import StringIO
import re
from typing import Any, Dict, List, Optional, Sequence, Tuple

from ase import Atoms as ASEAtoms
from ase.io import read as ase_read
from pymatgen.io.pwscf import PWInput

from app.schemas.common import Atom, Lattice, QeParameters, Structure, Vector3


def _ase_atoms_from_content(content: str) -> ASEAtoms:
    atoms_result = ase_read(StringIO(content), format="espresso-in")
    if isinstance(atoms_result, list):
        if not atoms_result:
            raise ValueError("ASEが構造を返しませんでした。")
        return atoms_result[0]
    return atoms_result


def _pymatgen_atoms_from_content(content: str) -> ASEAtoms:
    pw_input = PWInput.from_str(content)
    structure = pw_input.structure
    symbols = [
        getattr(site.specie, "symbol", str(site.specie))
        for site in structure.sites
    ]
    positions = [site.coords for site in structure.sites]
    return ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=structure.lattice.matrix.tolist(),
        pbc=True,
    )


def structure_from_ase(atoms_obj: ASEAtoms) -> Structure:
    symbols = atoms_obj.get_chemical_symbols()
    positions = atoms_obj.get_positions()
    parsed: List[Atom] = []
    for symbol, (x, y, z) in zip(symbols, positions):
        parsed.append(Atom(symbol=symbol, x=float(x), y=float(y), z=float(z)))
    lattice = _lattice_from_vectors(atoms_obj.get_cell().tolist())
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
    if isinstance(value, dict):
        return dict(value)
    try:
        return dict(value)  # type: ignore[call-overload,arg-type]
    except Exception:
        return {}


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
            return result if isinstance(result, dict) else None
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
