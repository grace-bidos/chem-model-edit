from __future__ import annotations

import pytest
from fastapi.testclient import TestClient

from main import app
from services.transplant import read_qe_input

CLIENT = TestClient(app)

SMALL_IN = """
&CONTROL
  calculation='relax'
/
&SYSTEM
  ibrav=0, nat=3, ntyp=2
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
 O 15.999 O.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
 H 1.0 0.0 0.0 0 0 0
 O 2.0 0.0 0.0 1 1 1
""".strip()

SMALL_OUT = """
...
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0
 H 1.0 0.0 0.0
 O 2.0 0.0 0.0
...
ATOMIC_POSITIONS angstrom
 H 0.1 0.0 0.0
 H 1.0 0.0 0.0
 O 2.0 -0.2 0.0
""".strip()

LARGE_IN = """
&CONTROL
  calculation='scf'
/
&SYSTEM
  ibrav=0, nat=4, ntyp=3
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus.UPF
 O 15.999 O.pbe-rrkjus.UPF
 C 12.0107 C.pbe-rrkjus.UPF
CELL_PARAMETERS angstrom
  5.0 0.0 0.0
  0.0 5.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
 H 1.0 0.0 0.0 0 0 0
 O 2.0 0.0 0.0 1 1 1
 C 3.0 0.0 0.0 0 0 0
""".strip()


def test_transplant_delta_ok():
    response = CLIENT.post(
        "/api/transforms/delta-transplant",
        json={"small_in": SMALL_IN, "small_out": SMALL_OUT, "large_in": LARGE_IN},
    )
    assert response.status_code == 200
    content = response.json()["content"]
    structure = read_qe_input(content)
    assert structure.atoms[0].x == pytest.approx(0.1)
    assert structure.atoms[0].y == pytest.approx(0.0)
    assert structure.atoms[2].y == pytest.approx(4.8)


def test_transplant_delta_missing_flags():
    small_in = SMALL_IN.replace("1 1 1", "").replace("0 0 0", "")
    response = CLIENT.post(
        "/api/transforms/delta-transplant",
        json={"small_in": small_in, "small_out": SMALL_OUT, "large_in": LARGE_IN},
    )
    assert response.status_code == 400
    message = response.json().get("error", {}).get("message", "")
    assert "フラグ" in message


def test_transplant_delta_missing_match():
    large_in = LARGE_IN.replace("2.0 0.0 0.0", "2.1 0.0 0.0", 1)
    response = CLIENT.post(
        "/api/transforms/delta-transplant",
        json={"small_in": SMALL_IN, "small_out": SMALL_OUT, "large_in": large_in},
    )
    assert response.status_code == 400
    message = response.json().get("error", {}).get("message", "")
    assert "一致する原子" in message


def test_transplant_delta_duplicate_match():
    large_in = LARGE_IN.replace(
        "C 3.0 0.0 0.0 0 0 0",
        "O 2.0 0.0 0.0 0 0 0",
    )
    response = CLIENT.post(
        "/api/transforms/delta-transplant",
        json={"small_in": SMALL_IN, "small_out": SMALL_OUT, "large_in": large_in},
    )
    assert response.status_code == 400
    message = response.json().get("error", {}).get("message", "")
    assert "複数" in message
