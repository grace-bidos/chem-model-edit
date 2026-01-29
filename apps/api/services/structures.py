from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from threading import RLock
from typing import Dict, Optional
from uuid import uuid4

from ase import Atoms as ASEAtoms

from models import QeParameters, Structure
from services.cif import atoms_to_cif, cif_to_bcif
from services.parse import extract_qe_params, parse_qe_atoms, structure_from_ase

@dataclass
class StoredStructure:
    atoms: ASEAtoms
    source: str
    cif: str
    params: QeParameters | None
    raw_input: str | None
    created_at: datetime


class StructureStore:
    def __init__(self) -> None:
        self._items: Dict[str, StoredStructure] = {}
        self._lock = RLock()

    def create(
        self,
        atoms: ASEAtoms,
        source: str,
        cif: str,
        params: QeParameters | None,
        raw_input: str | None,
    ) -> str:
        now = datetime.now(timezone.utc)
        structure_id = uuid4().hex
        entry = StoredStructure(
            atoms=atoms,
            source=source,
            cif=cif,
            params=params,
            raw_input=raw_input,
            created_at=now,
        )
        with self._lock:
            self._items[structure_id] = entry
        return structure_id

    def get(self, structure_id: str) -> Optional[StoredStructure]:
        with self._lock:
            entry = self._items.get(structure_id)
            if not entry:
                return None
            return entry


_STORE = StructureStore()


def create_structure_from_qe(
    content: str,
) -> tuple[str, Structure, str, QeParameters | None]:
    atoms, source = parse_qe_atoms(content)
    cif = atoms_to_cif(atoms)
    structure = structure_from_ase(atoms)
    params = extract_qe_params(content)
    structure_id = _STORE.create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=params,
        raw_input=content,
    )
    return structure_id, structure, source, params


def get_structure_bcif(
    structure_id: str,
    *,
    lossy: bool,
    precision: Optional[int],
) -> bytes:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return cif_to_bcif(entry.cif, lossy=lossy, precision=precision)


def get_structure_entry(structure_id: str) -> StoredStructure:
    entry = _STORE.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return entry


def get_structure(structure_id: str) -> Structure:
    entry = get_structure_entry(structure_id)
    return structure_from_ase(entry.atoms)


def get_structure_params(structure_id: str) -> QeParameters | None:
    entry = get_structure_entry(structure_id)
    return entry.params


def get_structure_raw_input(structure_id: str) -> str | None:
    entry = get_structure_entry(structure_id)
    return entry.raw_input


def register_structure_atoms(atoms: ASEAtoms, source: str) -> str:
    cif = atoms_to_cif(atoms)
    return _STORE.create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=None,
        raw_input=None,
    )
