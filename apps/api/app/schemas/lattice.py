from __future__ import annotations

from typing import Literal, Optional

from pydantic import Field, model_validator

from .base import ApiModel
from .common import Lattice, LatticeParams


class LatticeConvertRequest(ApiModel):
    from_: Literal["vectors", "params"] = Field(..., alias="from")
    unit: str = "angstrom"
    lattice: Optional[Lattice] = None
    params: Optional[LatticeParams] = None

    @model_validator(mode="after")
    def validate_source(self) -> "LatticeConvertRequest":
        if self.from_ == "vectors" and self.lattice is None:
            raise ValueError("lattice is required when from=vectors")
        if self.from_ == "params" and self.params is None:
            raise ValueError("params is required when from=params")
        return self


class LatticeConvertResponse(ApiModel):
    lattice: Lattice
    params: LatticeParams
    unit: str
