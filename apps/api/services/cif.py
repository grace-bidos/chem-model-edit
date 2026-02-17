from __future__ import annotations

from importlib import import_module
from io import BytesIO
from typing import Callable, cast

from ase import Atoms as ASEAtoms

ase_write = cast(Callable[..., None], getattr(import_module("ase.io"), "write"))


def atoms_to_cif(atoms: ASEAtoms, *, wrap: bool = False) -> str:
    buffer = BytesIO()
    ase_write(buffer, atoms, format="cif", wrap=wrap)
    return buffer.getvalue().decode("utf-8")
