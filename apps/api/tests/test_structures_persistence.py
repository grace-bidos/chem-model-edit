from __future__ import annotations

from pathlib import Path

from ase import Atoms as ASEAtoms
from fastapi.testclient import TestClient

import main
from services.cif import atoms_to_cif
from services import structures as structure_service
from services.structures import SQLiteStructureStore

QE_INPUT = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
/
&ELECTRONS
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 5.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0
 H 0.0 0.0 0.74
""".strip()


def test_sqlite_structure_store_persists_across_reopen(tmp_path: Path) -> None:
    db_path = tmp_path / "structures.sqlite3"
    atoms = ASEAtoms(symbols=["H", "H"], positions=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.74)])
    cif = atoms_to_cif(atoms)

    first = SQLiteStructureStore(db_path, reset_on_start=True)
    structure_id = first.create(
        atoms=atoms,
        source="qe",
        cif=cif,
        params=None,
        raw_input="raw-input",
    )

    reopened = SQLiteStructureStore(db_path, reset_on_start=False)
    entry = reopened.get(structure_id)

    assert entry is not None
    assert entry.source == "qe"
    assert entry.raw_input == "raw-input"
    assert len(entry.atoms) == 2


def test_structure_api_survives_store_reload(monkeypatch, tmp_path: Path) -> None:
    db_path = tmp_path / "structures.sqlite3"
    monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "sqlite")
    monkeypatch.setenv("STRUCTURE_STORE_DB_PATH", str(db_path))
    monkeypatch.setenv("STRUCTURE_STORE_RESET_ON_START", "1")
    structure_service.reload_structure_store()

    try:
        client = TestClient(main.app)
        created = client.post("/api/structures", json={"content": QE_INPUT})
        assert created.status_code == 200
        structure_id = created.json()["structure_id"]

        monkeypatch.setenv("STRUCTURE_STORE_RESET_ON_START", "0")
        structure_service.reload_structure_store()

        fetched = client.get(f"/api/structures/{structure_id}")
        assert fetched.status_code == 200
        payload = fetched.json()
        assert payload["raw_input"] == QE_INPUT
        assert len(payload["structure"]["atoms"]) == 2

        monkeypatch.setenv("STRUCTURE_STORE_RESET_ON_START", "1")
        structure_service.reload_structure_store()

        missing = client.get(f"/api/structures/{structure_id}")
        assert missing.status_code == 404
    finally:
        monkeypatch.setenv("STRUCTURE_STORE_BACKEND", "memory")
        monkeypatch.delenv("STRUCTURE_STORE_DB_PATH", raising=False)
        monkeypatch.delenv("STRUCTURE_STORE_RESET_ON_START", raising=False)
        structure_service.reload_structure_store()
