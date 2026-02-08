from __future__ import annotations

from typing import List, Literal, Optional

from pydantic import Field, model_validator

from .base import ApiModel
from .common import Structure


class SupercellGridAxis(ApiModel):
    row: Literal["a", "b"]
    col: Literal["a", "b"]

    @model_validator(mode="after")
    def check_distinct_axes(self) -> "SupercellGridAxis":
        if self.row == self.col:
            raise ValueError("row and col must map to different axes")
        return self


class SupercellGrid(ApiModel):
    rows: int = Field(..., ge=1)
    cols: int = Field(..., ge=1)
    tiles: List[List[str]]
    axis: Optional[SupercellGridAxis] = None

    @model_validator(mode="after")
    def check_tile_dimensions(self) -> "SupercellGrid":
        if len(self.tiles) != self.rows:
            raise ValueError(f"tiles has {len(self.tiles)} rows, expected {self.rows}")
        for index, row in enumerate(self.tiles):
            if len(row) != self.cols:
                raise ValueError(
                    f"tiles[{index}] has {len(row)} cols, expected {self.cols}"
                )
        return self


class SupercellBuildOptions(ApiModel):
    check_overlap: bool = False
    overlap_tolerance: Optional[float] = None
    validate_lattice: Literal["none", "warn", "error"] = "none"


class SupercellBuildOutput(ApiModel):
    include_structure: bool = False


class SupercellBuildRequest(ApiModel):
    base_structure_id: str
    grid: SupercellGrid
    options: Optional[SupercellBuildOptions] = None
    output: Optional[SupercellBuildOutput] = None


class SupercellBuildMeta(ApiModel):
    rows: int
    cols: int
    tile_count: int
    overlap_count: Optional[int] = None
    base_structure_id: str
    structure_ids_used: List[str]


class SupercellBuildResponse(ApiModel):
    structure_id: str
    structure: Optional[Structure] = None
    meta: SupercellBuildMeta
