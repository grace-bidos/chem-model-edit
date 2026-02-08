from __future__ import annotations

from pydantic import BaseModel, ConfigDict

class ApiModel(BaseModel):
    model_config = ConfigDict(
        populate_by_name=True,
        from_attributes=True,
        extra="forbid",
    )
