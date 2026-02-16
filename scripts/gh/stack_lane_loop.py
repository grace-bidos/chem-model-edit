#!/usr/bin/env python3
"""Lane-owner entrypoint for stack sync + PR readiness loop.

This script keeps responsibilities explicit:
- `gt` handles stack operations (`sync` before readiness checks).
- `gh` helper scripts handle readiness summary and optional watch/merge loop.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], dry_run: bool, allowed_exit_codes: set[int] | None = None) -> int:
    printable = " ".join(shlex.quote(part) for part in cmd)
    print(f"$ {printable}")
    if dry_run:
        return 0
    proc = subprocess.run(cmd, check=False)
    allowed = allowed_exit_codes or set()
    if proc.returncode != 0 and proc.returncode not in allowed:
        raise RuntimeError(f"command failed ({proc.returncode}): {printable}")
    return proc.returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run optional gt sync, then PR readiness summary and optional pr-autoloop watch/merge."
        )
    )
    parser.add_argument("pr", help="PR number or PR URL")
    parser.add_argument(
        "--gt-sync",
        action="store_true",
        help="Run `gt sync` before PR readiness commands",
    )
    parser.add_argument(
        "--watch",
        action="store_true",
        help="Pass through to pr-autoloop.py for polling mode",
    )
    parser.add_argument(
        "--merge-when-ready",
        action="store_true",
        help="Pass through to pr-autoloop.py to merge automatically when ready",
    )
    parser.add_argument(
        "--merge-method",
        choices=["merge", "squash", "rebase"],
        default="merge",
        help="Merge method used with --merge-when-ready (default: merge)",
    )
    parser.add_argument(
        "--delete-branch",
        action="store_true",
        help="Pass through to pr-autoloop.py",
    )
    parser.add_argument(
        "--resolve-outdated-threads",
        action="store_true",
        help="Pass through to pr-autoloop.py",
    )
    parser.add_argument(
        "--interval",
        type=int,
        default=20,
        help="Polling interval seconds for pr-autoloop.py (default: 20)",
    )
    parser.add_argument(
        "--max-wait",
        type=int,
        default=0,
        help="Max wait seconds in --watch mode (default: 0, no timeout)",
    )
    parser.add_argument(
        "--skip-summary",
        action="store_true",
        help="Skip one-shot readiness summary",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    pr_readiness = script_dir / "pr_readiness.py"
    pr_autoloop = script_dir / "pr-autoloop.py"

    try:
        autoloop_needed = (
            args.watch
            or args.merge_when_ready
            or args.delete_branch
            or args.resolve_outdated_threads
            or args.interval != 20
            or args.max_wait != 0
            or args.merge_method != "merge"
        )

        if args.gt_sync:
            run(["gt", "sync"], dry_run=args.dry_run)

        if not args.skip_summary:
            allowed = {2} if autoloop_needed else None
            run([str(pr_readiness), args.pr], dry_run=args.dry_run, allowed_exit_codes=allowed)

        if autoloop_needed:
            cmd = [str(pr_autoloop), args.pr]
            if args.watch:
                cmd.append("--watch")
            if args.merge_when_ready:
                cmd.append("--merge-when-ready")
            if args.delete_branch:
                cmd.append("--delete-branch")
            if args.resolve_outdated_threads:
                cmd.append("--resolve-outdated-threads")
            if args.interval != 20:
                cmd.extend(["--interval", str(args.interval)])
            if args.max_wait != 0:
                cmd.extend(["--max-wait", str(args.max_wait)])
            if args.merge_method != "merge":
                cmd.extend(["--merge-method", args.merge_method])
            run(cmd, dry_run=args.dry_run)
    except Exception as exc:  # noqa: BLE001
        print(f"error: {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
