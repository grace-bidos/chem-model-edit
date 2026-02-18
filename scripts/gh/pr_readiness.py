#!/usr/bin/env python3
"""Summarize PR readiness (checks, unresolved review threads, merge state)."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from typing import Any


def run_gh(args: list[str]) -> str:
    proc = subprocess.run(
        ["gh", *args],
        check=False,
        text=True,
        capture_output=True,
    )
    if proc.returncode != 0:
        cmd = "gh " + " ".join(shlex.quote(arg) for arg in args)
        raise RuntimeError(f"command failed ({proc.returncode}): {cmd}\n{proc.stderr.strip()}")
    return proc.stdout


def parse_pr_number(value: str) -> int:
    if value.isdigit():
        return int(value)
    if "/pull/" in value:
        tail = value.rsplit("/pull/", 1)[1].strip("/")
        if tail.isdigit():
            return int(tail)
    raise ValueError(f"invalid PR identifier: {value}")


def get_repo_owner_name() -> tuple[str, str]:
    out = run_gh(["repo", "view", "--json", "nameWithOwner"])
    payload = json.loads(out)
    owner, repo = payload["nameWithOwner"].split("/", 1)
    return owner, repo


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


def get_required_status_checks(base_ref: str) -> set[str] | None:
    owner, repo = get_repo_owner_name()
    try:
        out = run_gh(
            [
                "api",
                f"repos/{owner}/{repo}/branches/{base_ref}/protection",
                "--jq",
                ".required_status_checks.contexts[]?",
            ]
        )
    except Exception:
        return None

    contexts = {line.strip() for line in out.splitlines() if line.strip()}
    return contexts or None


def get_unresolved_thread_count(pr_number: int) -> int:
    owner, repo = get_repo_owner_name()
    query = """
query($owner:String!, $repo:String!, $num:Int!, $after:String) {
  repository(owner:$owner, name:$repo) {
    pullRequest(number:$num) {
      reviewThreads(first:100, after:$after) {
        nodes {
          isResolved
        }
        pageInfo {
          hasNextPage
          endCursor
        }
      }
    }
  }
}
""".strip()

    unresolved = 0
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
        payload = json.loads(out)
        review_threads = payload["data"]["repository"]["pullRequest"]["reviewThreads"]
        for node in review_threads.get("nodes", []):
            if not node.get("isResolved", False):
                unresolved += 1

        page_info = review_threads.get("pageInfo") or {}
        if not page_info.get("hasNextPage"):
            break
        cursor = page_info.get("endCursor")
        if not cursor:
            break

    return unresolved


def classify_checks(
    rollup: list[dict[str, Any]], required_checks: set[str] | None
) -> tuple[list[str], list[str]]:
    pending: list[str] = []
    failing: list[str] = []

    for item in rollup:
        item_type = item.get("__typename")
        name = item.get("name") or item.get("context") or "(unknown)"
        if required_checks is not None and name not in required_checks:
            continue
        if item_type == "CheckRun":
            status = item.get("status", "")
            conclusion = item.get("conclusion", "")
            if status != "COMPLETED":
                pending.append(name)
                continue
            if conclusion in {"SUCCESS", "NEUTRAL", "SKIPPED"}:
                continue
            failing.append(name)
            continue

        if item_type == "StatusContext":
            state = item.get("state", "")
            if state == "SUCCESS":
                continue
            if state in {"PENDING", "EXPECTED"}:
                pending.append(name)
                continue
            failing.append(name)

    return pending, failing


def is_ready(pr: dict[str, Any], unresolved_count: int, pending: list[str], failing: list[str]) -> bool:
    return (
        pr.get("state") == "OPEN"
        and not pr.get("isDraft", False)
        and pr.get("mergeStateStatus") == "CLEAN"
        and unresolved_count == 0
        and len(pending) == 0
        and len(failing) == 0
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize PR readiness for CI/review flow")
    parser.add_argument("pr", help="PR number or PR URL")
    args = parser.parse_args()

    try:
        pr_number = parse_pr_number(args.pr)
        pr = get_pr(pr_number)
        required_checks = get_required_status_checks(pr.get("baseRefName", "main"))
        unresolved_count = get_unresolved_thread_count(pr_number)
        pending_checks, failing_checks = classify_checks(
            pr.get("statusCheckRollup") or [],
            required_checks,
        )
    except Exception as exc:  # noqa: BLE001
        print(str(exc), file=sys.stderr)
        return 3

    print(f"PR #{pr['number']}: {pr['title']}")
    print(f"URL: {pr['url']}")
    print(
        "State:"
        f" {pr.get('state')}"
        f", Draft: {pr.get('isDraft', False)}"
        f", mergeStateStatus: {pr.get('mergeStateStatus')}"
        f", Branch: {pr.get('headRefName')} -> {pr.get('baseRefName')}"
    )
    print(
        "Checks:"
        f" pending={len(pending_checks)}"
        f", failing={len(failing_checks)}"
        f", unresolved_review_threads={unresolved_count}"
    )
    if required_checks is None:
        print("Required check scope: unavailable (fallback to all checks in rollup)")
    else:
        names = ", ".join(sorted(required_checks)) if required_checks else "(none)"
        print(f"Required check scope: {names}")

    if failing_checks:
        print("Failing checks:")
        for name in failing_checks:
            print(f"- {name}")

    if pending_checks:
        print("Pending checks:")
        for name in pending_checks:
            print(f"- {name}")

    merge_state = pr.get("mergeStateStatus")
    if merge_state == "BEHIND":
        print(
            "Action: branch is behind base. Rebase/merge base into head, then push and re-run checks."
        )

    if pr.get("state") == "MERGED":
        print("Result: already merged")
        return 0

    if pr.get("state") == "CLOSED":
        print("Result: closed without merge")
        return 1

    if is_ready(pr, unresolved_count, pending_checks, failing_checks):
        print("Result: merge-ready")
        return 0

    print("Result: not ready")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
