from __future__ import annotations

from typing import Any, Dict, List, Optional

from pydantic import Field

from .base import ApiModel


class Atom(ApiModel):
    symbol: str = Field(..., min_length=1, max_length=3)
    x: float
    y: float
    z: float


class Vector3(ApiModel):
    x: float
    y: float
    z: float


class Lattice(ApiModel):
    a: Vector3
    b: Vector3
    c: Vector3


class LatticeParams(ApiModel):
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


class Structure(ApiModel):
    atoms: List[Atom]
    lattice: Optional[Lattice] = None


class QeParameters(ApiModel):
    control: Dict[str, Any] = Field(default_factory=dict)
    system: Dict[str, Any] = Field(default_factory=dict)
    electrons: Dict[str, Any] = Field(default_factory=dict)
    ions: Dict[str, Any] = Field(default_factory=dict)
    cell: Dict[str, Any] = Field(default_factory=dict)
    pseudopotentials: Dict[str, str] = Field(default_factory=dict)
    kpoints: Optional[Dict[str, Any]] = None


class Pagination(ApiModel):
    total: int
    limit: int
    offset: int
