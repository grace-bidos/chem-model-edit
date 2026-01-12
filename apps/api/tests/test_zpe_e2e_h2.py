from __future__ import annotations

import math
import os
from pathlib import Path

import fakeredis
import pytest

from services.zpe import settings as zpe_settings
from services.zpe import worker
from services.zpe.result_store import RedisResultStore


QE_INPUT_TEMPLATE = """
&CONTROL
  calculation='scf'
  prefix='ase_calc'
  outdir='calc_scratch'
/
&SYSTEM
  ibrav=0, nat=2, ntyp=1
  ecutwfc=10.0, ecutrho=80.0
  occupations='smearing'
  smearing='mp'
  degauss=0.02
/
&ELECTRONS
  conv_thr=1.0d-4
  electron_maxstep=30
/
ATOMIC_SPECIES
 H 1.0079 {pseudo}
CELL_PARAMETERS angstrom
  6.0 0.0 0.0
  0.0 6.0 0.0
  0.0 0.0 6.0
ATOMIC_POSITIONS angstrom
 H 0.0 0.0 0.0 1 1 1
 H 0.0 0.0 0.74 1 1 1
K_POINTS gamma
""".strip()


@pytest.mark.e2e
def test_zpe_e2e_h2(monkeypatch):
    if os.environ.get("ZPE_E2E") != "1":
        pytest.skip("Set ZPE_E2E=1 to run QE E2E test")

    pw_path = os.environ.get("ZPE_PW_PATH")
    pseudo_dir = os.environ.get("ZPE_PSEUDO_DIR")
    pseudo_name = os.environ.get("ZPE_E2E_PSEUDO", "H.pbe-rrkjus_psl.1.0.0.UPF")

    if not pw_path or not Path(pw_path).exists():
        pytest.skip("ZPE_PW_PATH is missing or invalid")
    if not pseudo_dir or not Path(pseudo_dir).is_dir():
        pytest.skip("ZPE_PSEUDO_DIR is missing or invalid")
    if not (Path(pseudo_dir) / pseudo_name).exists():
        pytest.skip("Pseudo file not found in ZPE_PSEUDO_DIR")

    monkeypatch.setenv("ZPE_PW_PATH", pw_path)
    monkeypatch.setenv("ZPE_PSEUDO_DIR", pseudo_dir)
    monkeypatch.setenv("ZPE_USE_MPI", "false")
    monkeypatch.setenv("ZPE_NP_CORE", "1")
    zpe_settings.get_zpe_settings.cache_clear()

    fake = fakeredis.FakeRedis()
    monkeypatch.setattr(worker, "get_result_store", lambda: RedisResultStore(redis=fake))

    qe_input = QE_INPUT_TEMPLATE.format(pseudo=pseudo_name)
    payload = {
        "content": qe_input,
        "mobile_indices": [0, 1],
        "use_environ": False,
        "input_dir": None,
        "calc_mode": "new",
        "job_id": "e2e-h2",
    }

    result = worker.run_zpe_job(payload)
    freqs = result["freqs_cm"]
    assert len(freqs) == 6
    assert all(math.isfinite(freq) for freq in freqs)

    low_cut = float(result["low_cut_cm"])
    positive = [freq for freq in freqs if math.isfinite(freq) and freq > low_cut]
    max_freq = max(positive)
    ref_freq = 4401.0
    tol_freq = float(os.environ.get("ZPE_E2E_FREQ_TOL_FRAC", "0.25"))
    assert abs(max_freq - ref_freq) <= ref_freq * tol_freq

    ref_zpe = 0.26
    tol_zpe = float(os.environ.get("ZPE_E2E_ZPE_TOL_FRAC", "0.35"))
    zpe_ev = result["zpe_ev"]
    assert abs(zpe_ev - ref_zpe) <= ref_zpe * tol_zpe
