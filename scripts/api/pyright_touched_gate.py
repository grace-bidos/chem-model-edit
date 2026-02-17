#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any


def _run(cmd: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(cmd, text=True, cwd=cwd)


def _to_repo_relative(path_text: str, repo_root: Path) -> str | None:
    raw = Path(path_text)
    if not raw.is_absolute():
        raw = (repo_root / raw).resolve()
    try:
        return raw.resolve().relative_to(repo_root).as_posix()
    except ValueError:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-ref", default="origin/main")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[2]
    api_root = repo_root / "apps" / "api"

    diff_names = _run(
        [
            "git",
            "diff",
            "--name-only",
            f"{args.base_ref}...HEAD",
            "--",
            "apps/api",
        ],
        cwd=repo_root,
    )
    touched_files = {
        line.strip()
        for line in diff_names.splitlines()
        if line.strip().endswith(".py")
    }

    if not touched_files:
        print(
            f"Pyright touched-files gate: no changed Python files under apps/api vs {args.base_ref}; gate passes."
        )
        return 0

    try:
        raw_json = _run(
            ["uv", "run", "pyright", "--project", "pyrightconfig.json", "--outputjson"],
            cwd=api_root,
        )
    except subprocess.CalledProcessError as exc:
        output = exc.stdout if isinstance(exc.stdout, str) else ""
        if not output:
            print(
                "Pyright touched-files gate: failed to capture pyright JSON output.",
                file=sys.stderr,
            )
            return 2
        raw_json = output
    payload: dict[str, Any] = json.loads(raw_json)
    diagnostics: list[dict[str, Any]] = payload.get("generalDiagnostics", [])

    warning_map: dict[str, list[tuple[int, str, str]]] = defaultdict(list)
    for diag in diagnostics:
        if diag.get("severity") != "warning":
            continue
        file_text = diag.get("file")
        if not isinstance(file_text, str):
            continue
        repo_rel = _to_repo_relative(file_text, repo_root)
        if repo_rel is None or repo_rel not in touched_files:
            continue

        range_info = diag.get("range", {})
        start = range_info.get("start", {})
        line = int(start.get("line", 0)) + 1
        rule = str(diag.get("rule") or "unknown-rule")
        message = str(diag.get("message") or "").strip()
        warning_map[repo_rel].append((line, rule, message))

    warning_total = sum(len(items) for items in warning_map.values())
    if warning_total == 0:
        print(
            f"Pyright touched-files gate: 0 warnings across {len(touched_files)} changed Python files (base {args.base_ref})."
        )
        return 0

    print(
        f"Pyright touched-files gate FAILED: {warning_total} warnings in {len(warning_map)} changed files (base {args.base_ref})."
    )
    for filename in sorted(warning_map):
        print(f"- {filename}")
        for line, rule, message in sorted(warning_map[filename]):
            print(f"  - L{line} [{rule}] {message}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
