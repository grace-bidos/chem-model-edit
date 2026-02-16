from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from app.api import app


def _render_schema() -> str:
    schema = app.openapi()
    return json.dumps(schema, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export FastAPI OpenAPI schema to api-client artifacts.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if the checked-in OpenAPI artifact is out of date.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[3]
    output_path = repo_root / "packages" / "api-client" / "openapi" / "openapi.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rendered = _render_schema()

    if args.check:
        if not output_path.exists():
            print(f"missing {output_path}")
            return 1
        existing = output_path.read_text(encoding="utf-8")
        if existing != rendered:
            print(
                "OpenAPI artifact drift detected. Run:\n"
                "  PYTHONPATH=apps/api uv run --project apps/api python "
                "apps/api/scripts/export_openapi.py"
            )
            return 1
        print(f"up-to-date {output_path}")
        return 0

    output_path.write_text(rendered, encoding="utf-8")
    print(f"wrote {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
