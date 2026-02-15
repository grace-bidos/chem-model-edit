#!/usr/bin/env python3
"""Watch a GitHub PR and optionally merge it when it becomes ready.

Requirements:
- gh CLI authenticated (`gh auth status`)
- repository context available from current working directory
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any


@dataclass
class ThreadInfo:
    thread_id: str
    is_outdated: bool
    author: str
    url: str


def run_gh(args: list[str]) -> str:
    proc = subprocess.run(
        ["gh", *args],
        check=False,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        cmd = "gh " + " ".join(shlex.quote(a) for a in args)
        raise RuntimeError(f"command failed ({proc.returncode}): {cmd}\n{proc.stderr.strip()}")
    return proc.stdout


def get_repo_owner_name() -> tuple[str, str]:
    out = run_gh(["repo", "view", "--json", "nameWithOwner"]) 
    data = json.loads(out)
    owner, repo = data["nameWithOwner"].split("/", 1)
    return owner, repo


def parse_pr_number(pr: str) -> int:
    if pr.isdigit():
        return int(pr)
    if "/pull/" in pr:
        tail = pr.rsplit("/pull/", 1)[1].strip("/")
        if tail.isdigit():
            return int(tail)
    raise ValueError(f"invalid PR identifier: {pr}")


def get_pr(pr_number: int) -> dict[str, Any]:
    fields = [
        "number",
        "title",
        "url",
        "state",
        "isDraft",
        "mergeStateStatus",
        "headRefName",
        "baseRefName",
        "statusCheckRollup",
    ]
    out = run_gh(["pr", "view", str(pr_number), "--json", ",".join(fields)])
    return json.loads(out)


def get_unresolved_threads(pr_number: int) -> list[ThreadInfo]:
    owner, repo = get_repo_owner_name()
    query = (
        "query($owner:String!,$repo:String!,$num:Int!,$after:String){"
        "repository(owner:$owner,name:$repo){"
        "pullRequest(number:$num){"
        "reviewThreads(first:100, after:$after){"
        "nodes{"
        "id isResolved isOutdated "
        "comments(first:1){nodes{author{login} url}}"
        "} pageInfo{hasNextPage endCursor}"
        "}}}}}"
    )
    unresolved: list[ThreadInfo] = []
    cursor: str | None = None

    while True:
        args = [
            "api",
            "graphql",
            "-f",
            f"query={query}",
            "-F",
            f"owner={owner}",
            "-F",
            f"repo={repo}",
            "-F",
            f"num={pr_number}",
        ]
        if cursor:
            args.extend(["-F", f"after={cursor}"])
        else:
            args.extend(["-F", "after="])

        out = run_gh(args)
        data = json.loads(out)
        threads = data["data"]["repository"]["pullRequest"]["reviewThreads"]
        nodes = threads.get("nodes", [])

        for node in nodes:
            if node["isResolved"]:
                continue
            comment_nodes = node.get("comments", {}).get("nodes", [])
            if not comment_nodes:
                continue
            first = comment_nodes[0]
            unresolved.append(
                ThreadInfo(
                    thread_id=node["id"],
                    is_outdated=bool(node.get("isOutdated", False)),
                    author=(first.get("author") or {}).get("login", "unknown"),
                    url=first.get("url", ""),
                )
            )

        page_info = threads.get("pageInfo") or {}
        if not page_info.get("hasNextPage"):
            break
        cursor = page_info.get("endCursor")
        if not cursor:
            break

    return unresolved


def classify_checks(rollup: list[dict[str, Any]]) -> tuple[list[str], list[str]]:
    pending: list[str] = []
    failing: list[str] = []

    for item in rollup:
        name = item.get("name") or item.get("context") or "(unknown)"
        kind = item.get("__typename")
        if kind == "CheckRun":
            status = item.get("status", "")
            conclusion = item.get("conclusion", "")
            if status != "COMPLETED":
                pending.append(name)
                continue
            if conclusion in {"SUCCESS", "NEUTRAL", "SKIPPED"}:
                continue
            failing.append(name)
        elif kind == "StatusContext":
            state = item.get("state", "")
            if state in {"SUCCESS"}:
                continue
            if state in {"PENDING", "EXPECTED"}:
                pending.append(name)
            else:
                failing.append(name)

    return pending, failing


def can_attempt_merge(pr: dict[str, Any], unresolved_count: int) -> bool:
    if pr["state"] != "OPEN":
        return False
    if pr.get("isDraft"):
        return False
    if unresolved_count > 0:
        return False
    return pr.get("mergeStateStatus") == "CLEAN"


def merge_pr(pr_number: int, method: str, delete_branch: bool) -> None:
    args = ["pr", "merge", str(pr_number), f"--{method}"]
    if delete_branch:
        args.append("--delete-branch")
    run_gh(args)


def now_str() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%SZ")


def main() -> int:
    # Keep logs visible during long-running watch mode even when stdout is piped.
    try:
        sys.stdout.reconfigure(line_buffering=True)
        sys.stderr.reconfigure(line_buffering=True)
    except Exception:
        pass

    parser = argparse.ArgumentParser(
        description="Watch PR status and optionally auto-merge when ready."
    )
    parser.add_argument("pr", help="PR number or URL")
    parser.add_argument("--watch", action="store_true", help="Poll until terminal state")
    parser.add_argument("--interval", type=int, default=20, help="Polling interval seconds")
    parser.add_argument(
        "--max-wait",
        type=int,
        default=0,
        help="Max wait seconds in --watch mode (0 means no timeout)",
    )
    parser.add_argument(
        "--merge-when-ready",
        action="store_true",
        help="Merge automatically when PR becomes ready",
    )
    parser.add_argument(
        "--merge-method",
        choices=["merge", "squash", "rebase"],
        default="merge",
        help="Merge method used with --merge-when-ready",
    )
    parser.add_argument(
        "--delete-branch",
        action="store_true",
        help="Delete remote/local branch after merge",
    )
    parser.add_argument(
        "--resolve-outdated-threads",
        action="store_true",
        help="Resolve only outdated unresolved review threads automatically",
    )

    args = parser.parse_args()
    pr_number = parse_pr_number(args.pr)
    start = time.time()

    while True:
        pr = get_pr(pr_number)
        unresolved = get_unresolved_threads(pr_number)

        if args.resolve_outdated_threads:
            for thread in unresolved:
                if not thread.is_outdated:
                    continue
                run_gh(
                    [
                        "api",
                        "graphql",
                        "-f",
                        "query=mutation($id:ID!){resolveReviewThread(input:{threadId:$id}){thread{id}}}",
                        "-F",
                        f"id={thread.thread_id}",
                    ]
                )
            unresolved = get_unresolved_threads(pr_number)

        pending_checks, failing_checks = classify_checks(pr.get("statusCheckRollup") or [])

        print(f"[{now_str()}] PR #{pr['number']} {pr['title']}")
        print(f"  url: {pr['url']}")
        print(
            "  state:"
            f" {pr['state']}"
            f", draft={pr.get('isDraft', False)}"
            f", mergeStateStatus={pr.get('mergeStateStatus')}"
        )
        print(
            f"  checks: pending={len(pending_checks)}, failing={len(failing_checks)},"
            f" unresolved_threads={len(unresolved)}"
        )
        if failing_checks:
            print("  failing checks:")
            for name in failing_checks:
                print(f"    - {name}")
        if unresolved:
            print("  unresolved threads:")
            for t in unresolved:
                flag = " [outdated]" if t.is_outdated else ""
                print(f"    - {t.author}{flag}: {t.url}")

        if pr["state"] == "MERGED":
            print("PR is already merged.")
            return 0
        if pr["state"] == "CLOSED":
            print("PR is closed without merge.")
            return 1

        if args.merge_when_ready and can_attempt_merge(pr, len(unresolved)):
            print(
                f"Ready to merge (method={args.merge_method}, delete_branch={args.delete_branch})."
            )
            merge_pr(pr_number, args.merge_method, args.delete_branch)
            print("Merge command submitted.")
            return 0

        if not args.watch:
            if can_attempt_merge(pr, len(unresolved)):
                print("PR is merge-ready.")
                return 0
            return 2

        if args.max_wait > 0 and (time.time() - start) >= args.max_wait:
            print("Timed out waiting for PR readiness.", file=sys.stderr)
            return 3

        time.sleep(max(5, args.interval))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:  # pragma: no cover
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
