from __future__ import annotations

import math

from app.schemas.common import Lattice, LatticeParams, Vector3
from services.lattice import params_to_vectors, vectors_to_params


def test_vectors_to_params_orthogonal():
    lattice = Lattice(
        a=Vector3(x=5.0, y=0.0, z=0.0),
        b=Vector3(x=0.0, y=6.0, z=0.0),
        c=Vector3(x=0.0, y=0.0, z=7.0),
    )
    params = vectors_to_params(lattice)
    assert math.isclose(params.a, 5.0, rel_tol=1e-6)
    assert math.isclose(params.b, 6.0, rel_tol=1e-6)
    assert math.isclose(params.c, 7.0, rel_tol=1e-6)
    assert math.isclose(params.alpha, 90.0, rel_tol=1e-6)
    assert math.isclose(params.beta, 90.0, rel_tol=1e-6)
    assert math.isclose(params.gamma, 90.0, rel_tol=1e-6)


def test_params_to_vectors_roundtrip():
    params = LatticeParams(a=4.0, b=5.0, c=6.0, alpha=90.0, beta=90.0, gamma=120.0)
    lattice = params_to_vectors(params)
    roundtrip = vectors_to_params(lattice)
    assert math.isclose(roundtrip.a, params.a, rel_tol=1e-6)
    assert math.isclose(roundtrip.b, params.b, rel_tol=1e-6)
    assert math.isclose(roundtrip.c, params.c, rel_tol=1e-6)
    assert math.isclose(roundtrip.alpha, params.alpha, rel_tol=1e-6)
    assert math.isclose(roundtrip.beta, params.beta, rel_tol=1e-6)
    assert math.isclose(roundtrip.gamma, params.gamma, rel_tol=1e-6)
