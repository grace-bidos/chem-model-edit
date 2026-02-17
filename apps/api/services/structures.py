from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import import_module
from io import StringIO
import json
import os
from pathlib import Path
import sqlite3
from threading import RLock
from typing import Any, Callable, Dict, Optional, cast
from uuid import uuid4

from ase import Atoms as ASEAtoms
from app.schemas.common import QeParameters, Structure
from services.cif import atoms_to_cif
from services.parse import extract_qe_params, parse_qe_atoms, structure_from_ase

ase_read = cast(Callable[..., Any], getattr(import_module("ase.io"), "read"))


@dataclass
class StoredStructure:
    atoms: ASEAtoms
    source: str
    cif: str
    params: QeParameters | None
    raw_input: str | None
    created_at: datetime


class InMemoryStructureStore:
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
            return self._items.get(structure_id)


class SQLiteStructureStore:
    def __init__(self, db_path: str | Path, *, reset_on_start: bool = False) -> None:
        self._db_path = Path(db_path)
        self._lock = RLock()
        self._db_path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize_schema(reset_on_start=reset_on_start)

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(str(self._db_path))
        connection.row_factory = sqlite3.Row
        return connection

    def _initialize_schema(self, *, reset_on_start: bool) -> None:
        with self._lock, self._connect() as connection:
            if reset_on_start:
                connection.execute("DROP TABLE IF EXISTS structures")
            connection.execute(
                """
                CREATE TABLE IF NOT EXISTS structures (
                  structure_id TEXT PRIMARY KEY,
                  source TEXT NOT NULL,
                  cif TEXT NOT NULL,
                  params_json TEXT,
                  raw_input TEXT,
                  created_at TEXT NOT NULL
                )
                """
            )
            connection.commit()

    def create(
        self,
        atoms: ASEAtoms,
        source: str,
        cif: str,
        params: QeParameters | None,
        raw_input: str | None,
    ) -> str:
        structure_id = uuid4().hex
        now = datetime.now(timezone.utc).isoformat()
        params_json = params.model_dump_json() if params is not None else None
        with self._lock, self._connect() as connection:
            connection.execute(
                """
                INSERT INTO structures (
                  structure_id,
                  source,
                  cif,
                  params_json,
                  raw_input,
                  created_at
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (structure_id, source, cif, params_json, raw_input, now),
            )
            connection.commit()
        return structure_id

    def get(self, structure_id: str) -> Optional[StoredStructure]:
        with self._lock, self._connect() as connection:
            row = connection.execute(
                """
                SELECT structure_id, source, cif, params_json, raw_input, created_at
                FROM structures
                WHERE structure_id = ?
                """,
                (structure_id,),
            ).fetchone()
        if row is None:
            return None

        cif = cast(str, row["cif"])
        params_json = cast(Optional[str], row["params_json"])
        params: QeParameters | None = None
        if params_json:
            params = QeParameters.model_validate(json.loads(params_json))

        created_at_raw = cast(str, row["created_at"])
        created_at = datetime.fromisoformat(created_at_raw)
        if created_at.tzinfo is None:
            created_at = created_at.replace(tzinfo=timezone.utc)

        return StoredStructure(
            atoms=_atoms_from_cif(cif),
            source=cast(str, row["source"]),
            cif=cif,
            params=params,
            raw_input=cast(Optional[str], row["raw_input"]),
            created_at=created_at,
        )


def _atoms_from_cif(cif_text: str) -> ASEAtoms:
    atoms = ase_read(StringIO(cif_text), format="cif")
    if not isinstance(atoms, ASEAtoms):
        raise ValueError("failed to parse CIF into ASE atoms")
    return atoms


def _default_db_path() -> Path:
    cwd = Path.cwd()
    for base in (cwd, *cwd.parents):
        if (base / ".git").exists():
            return base / ".just-runtime" / "structures.sqlite3"
    return cwd / ".just-runtime" / "structures.sqlite3"


def _as_bool(value: str | None, default: bool = False) -> bool:
    if value is None:
        return default
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    return default


def _build_store() -> InMemoryStructureStore | SQLiteStructureStore:
    backend = os.getenv("STRUCTURE_STORE_BACKEND", "sqlite").strip().lower()
    if backend == "memory":
        return InMemoryStructureStore()
    if backend != "sqlite":
        raise ValueError("STRUCTURE_STORE_BACKEND must be 'sqlite' or 'memory'")

    db_path_raw = os.getenv("STRUCTURE_STORE_DB_PATH")
    db_path = Path(db_path_raw).expanduser() if db_path_raw else _default_db_path()
    reset_on_start = _as_bool(os.getenv("STRUCTURE_STORE_RESET_ON_START"), False)
    return SQLiteStructureStore(db_path, reset_on_start=reset_on_start)


_store = _build_store()


def reload_structure_store() -> None:
    global _store
    _store = _build_store()


def create_structure_from_qe(
    content: str,
) -> tuple[str, Structure, str, QeParameters | None]:
    atoms, source = parse_qe_atoms(content)
    cif = atoms_to_cif(atoms)
    structure = structure_from_ase(atoms)
    params = extract_qe_params(content)
    structure_id = _store.create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=params,
        raw_input=content,
    )
    return structure_id, structure, source, params


def get_structure_cif(structure_id: str) -> str:
    entry = _store.get(structure_id)
    if not entry:
        raise KeyError(structure_id)
    return entry.cif


def get_structure_entry(structure_id: str) -> StoredStructure:
    entry = _store.get(structure_id)
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
    return _store.create(
        atoms=atoms,
        source=source,
        cif=cif,
        params=None,
        raw_input=None,
    )
