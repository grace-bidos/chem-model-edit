from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, cast
from uuid import uuid4

import fakeredis
import pytest

from services.zpe import settings as zpe_settings
from services.zpe.result_store import RedisResultStore
from services.zpe.thermo import HC_EV_CM
from services.zpe import worker


@dataclass(frozen=True)
class MoleculeCase:
    name: str
    natoms: int
    qe_input: str
    targets_cm: Sequence[float]
    pseudos: Sequence[str]


def _ensure_env() -> tuple[str, str]:
    if os.environ.get("ZPE_E2E") != "1":
        pytest.skip("Set ZPE_E2E=1 to run QE E2E tests")

    pw_path = os.environ.get("ZPE_PW_PATH")
    pseudo_dir = os.environ.get("ZPE_PSEUDO_DIR")
    if not pw_path or not Path(pw_path).exists():
        pytest.skip("ZPE_PW_PATH is missing or invalid")
    if not pseudo_dir or not Path(pseudo_dir).is_dir():
        pytest.skip("ZPE_PSEUDO_DIR is missing or invalid")
    return pw_path, pseudo_dir


def _set_env(monkeypatch, pw_path: str, pseudo_dir: str) -> None:
    monkeypatch.setenv("ZPE_PW_PATH", pw_path)
    monkeypatch.setenv("ZPE_PSEUDO_DIR", pseudo_dir)
    monkeypatch.setenv("ZPE_USE_MPI", "false")
    monkeypatch.setenv("ZPE_NP_CORE", "1")
    zpe_settings.get_zpe_settings.cache_clear()


def _run_case(monkeypatch, case: MoleculeCase) -> dict[str, object]:
    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(worker, "get_result_store", lambda: RedisResultStore(redis=fake))

    job_id = f"e2e-{case.name}-{uuid4().hex[:8]}"
    payload = {
        "content": case.qe_input,
        "mobile_indices": list(range(case.natoms)),
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "new",
        "job_id": job_id,
    }
    return worker.run_zpe_job(payload)


def _filter_positive(freqs: Iterable[float], *, low_cut: float) -> List[float]:
    return [freq for freq in freqs if math.isfinite(freq) and freq > low_cut]


def _assert_targets(freqs: Sequence[float], targets: Sequence[float], tol_frac: float) -> None:
    for target in targets:
        closest = min(abs(freq - target) for freq in freqs)
        assert closest <= target * tol_frac


def _assert_zpe(result: dict[str, object], targets: Sequence[float], tol_frac: float) -> None:
    expected = 0.5 * HC_EV_CM * float(sum(targets))
    zpe_ev = float(cast(float, result["zpe_ev"]))
    assert abs(zpe_ev - expected) <= expected * tol_frac


CASES: Sequence[MoleculeCase] = [
    MoleculeCase(
        name="n2",
        natoms=2,
        qe_input=(
            """
&CONTROL
  calculation='scf'
  prefix='ase_calc'
  outdir='calc_scratch'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
  ecutwfc=30.0, ecutrho=240.0
  occupations='fixed'
/
&ELECTRONS
  conv_thr=1.0d-6
  electron_maxstep=80
/
ATOMIC_SPECIES
 N 14.0067 N.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  10.0 0.0 0.0
  0.0 10.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 N 0.0 0.0 0.0 1 1 1
 N 0.0 0.0 1.10 1 1 1
K_POINTS gamma
"""
        ).strip(),
        targets_cm=[2330.0],
        pseudos=["N.pbe-n-rrkjus_psl.1.0.0.UPF"],
    ),
    MoleculeCase(
        name="co",
        natoms=2,
        qe_input=(
            """
&CONTROL
  calculation='scf'
  prefix='ase_calc'
  outdir='calc_scratch'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=2
  ecutwfc=30.0, ecutrho=240.0
  occupations='fixed'
/
&ELECTRONS
  conv_thr=1.0d-6
  electron_maxstep=80
/
ATOMIC_SPECIES
 C 12.0107 C.pbe-n-rrkjus_psl.1.0.0.UPF
 O 15.999 O.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  10.0 0.0 0.0
  0.0 10.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 C 0.0 0.0 0.0 1 1 1
 O 0.0 0.0 1.13 1 1 1
K_POINTS gamma
"""
        ).strip(),
        targets_cm=[2143.0],
        pseudos=["C.pbe-n-rrkjus_psl.1.0.0.UPF", "O.pbe-n-rrkjus_psl.1.0.0.UPF"],
    ),
    MoleculeCase(
        name="co2",
        natoms=3,
        qe_input=(
            """
&CONTROL
  calculation='scf'
  prefix='ase_calc'
  outdir='calc_scratch'
/
&SYSTEM
  ibrav=0, nat=3, ntyp=2
  ecutwfc=30.0, ecutrho=240.0
  occupations='fixed'
/
&ELECTRONS
  conv_thr=1.0d-6
  electron_maxstep=80
/
ATOMIC_SPECIES
 C 12.0107 C.pbe-n-rrkjus_psl.1.0.0.UPF
 O 15.999 O.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  10.0 0.0 0.0
  0.0 10.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 O 0.0 0.0 -1.16 1 1 1
 C 0.0 0.0 0.0 1 1 1
 O 0.0 0.0 1.16 1 1 1
K_POINTS gamma
"""
        ).strip(),
        targets_cm=[667.0, 1388.0, 2349.0],
        pseudos=["C.pbe-n-rrkjus_psl.1.0.0.UPF", "O.pbe-n-rrkjus_psl.1.0.0.UPF"],
    ),
    MoleculeCase(
        name="h2o",
        natoms=3,
        qe_input=(
            """
&CONTROL
  calculation='scf'
  prefix='ase_calc'
  outdir='calc_scratch'
/
&SYSTEM
  ibrav=0, nat=3, ntyp=2
  ecutwfc=30.0, ecutrho=240.0
  occupations='fixed'
/
&ELECTRONS
  conv_thr=1.0d-6
  electron_maxstep=80
/
ATOMIC_SPECIES
 H 1.0079 H.pbe-rrkjus_psl.1.0.0.UPF
 O 15.999 O.pbe-n-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  10.0 0.0 0.0
  0.0 10.0 0.0
  0.0 0.0 10.0
ATOMIC_POSITIONS angstrom
 O 0.0000 0.0000 0.0000 1 1 1
 H 0.9572 0.0000 0.0000 1 1 1
 H -0.2390 0.9270 0.0000 1 1 1
K_POINTS gamma
"""
        ).strip(),
        targets_cm=[1595.0, 3657.0, 3756.0],
        pseudos=["H.pbe-rrkjus_psl.1.0.0.UPF", "O.pbe-n-rrkjus_psl.1.0.0.UPF"],
    ),
]


@pytest.mark.e2e
@pytest.mark.parametrize("case", CASES, ids=[case.name for case in CASES])
def test_zpe_e2e_molecules(monkeypatch, case: MoleculeCase):
    pw_path, pseudo_dir = _ensure_env()
    _set_env(monkeypatch, pw_path, pseudo_dir)

    for pseudo in case.pseudos:
        if not (Path(pseudo_dir) / pseudo).exists():
            pytest.skip(f"Pseudo file not found: {pseudo}")

    result = _run_case(monkeypatch, case)
    freqs = cast(List[float], result["freqs_cm"])
    assert len(freqs) == 3 * case.natoms

    low_cut = float(cast(float, result["low_cut_cm"]))
    positive = _filter_positive(freqs, low_cut=low_cut)

    tol_freq = float(os.environ.get("ZPE_E2E_FREQ_TOL_FRAC", "0.3"))
    _assert_targets(positive, case.targets_cm, tol_freq)

    tol_zpe = float(os.environ.get("ZPE_E2E_ZPE_TOL_FRAC", "0.5"))
    _assert_zpe(result, case.targets_cm, tol_zpe)
