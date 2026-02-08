from __future__ import annotations

import json
from pathlib import Path

from app.api import app


def main() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    output_path = repo_root / "packages" / "api-client" / "openapi" / "openapi.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    schema = app.openapi()
    output_path.write_text(
        json.dumps(schema, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
