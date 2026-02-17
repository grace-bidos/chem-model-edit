#!/usr/bin/env python3
"""Generate phase-1 security gate promotion evidence from GitHub Actions runs.

The script collects recent CI runs on the target branch, inspects API security steps,
parses uploaded security summary artifacts, and emits:
- Markdown report for humans
- JSON report for machine processing

Exit behavior:
- Default: exit 0 for policy findings (non-ready/failed baselines), non-zero only on script failures.
- Strict mode: exit non-zero when promotion criteria are not met.
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

GH_TIMEOUT_SECONDS = 60


@dataclass
class StepStatus:
    name: str
    status: str
    conclusion: str


@dataclass
class RunEvidence:
    run_id: int
    run_number: int
    created_at: str
    url: str
    status: str
    conclusion: str
    branch: str
    event: str
    api_job_conclusion: str | None
    bandit_step: str
    pip_audit_step: str
    bandit_high: int | None
    bandit_medium: int | None
    bandit_low: int | None
    bandit_total: int | None
    pip_audit_vulnerable_packages: int | None
    pip_audit_vulnerabilities: int | None
    baseline_pass: bool | None
    artifact_found: bool
    notes: list[str]


# Keep these patterns in sync with CI summary lines in `.github/workflows/ci.yml`:
# - "- Bandit: high={high}, medium={medium}, low={low}, total={len(issues)}"
# - "- Pip-audit: vulnerable_packages={pkg_count}, vulnerabilities={vuln_count}"
BANDIT_RE = re.compile(
    r"- Bandit: high=(?P<high>\d+), medium=(?P<medium>\d+), low=(?P<low>\d+), total=(?P<total>\d+)"
)
PIP_AUDIT_RE = re.compile(
    r"- Pip-audit: vulnerable_packages=(?P<pkg>\d+), vulnerabilities=(?P<vuln>\d+)"
)


class GhError(RuntimeError):
    """Raised when gh command fails."""


def run_gh(args: list[str]) -> str:
    cmd = "gh " + " ".join(shlex.quote(a) for a in args)
    try:
        proc = subprocess.run(
            ["gh", *args],
            check=False,
            text=True,
            capture_output=True,
            timeout=GH_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        partial_stdout = (exc.stdout or "").strip()
        partial_stderr = (exc.stderr or "").strip()
        details = "\n".join(part for part in [partial_stdout, partial_stderr] if part) or "(no output)"
        raise GhError(
            f"command timed out after {GH_TIMEOUT_SECONDS}s: {cmd}\n{details}"
        ) from exc

    if proc.returncode != 0:
        stderr = proc.stderr.strip() or "(no stderr)"
        raise GhError(f"command failed ({proc.returncode}): {cmd}\n{stderr}")
    return proc.stdout


def run_gh_json(args: list[str]) -> Any:
    out = run_gh(args)
    try:
        return json.loads(out)
    except json.JSONDecodeError as exc:
        cmd = "gh " + " ".join(shlex.quote(a) for a in args)
        raise GhError(f"invalid JSON from: {cmd}\n{exc}") from exc


def normalize(value: str | None) -> str:
    return (value or "unknown").strip().lower()


def step_map(steps: list[dict[str, Any]]) -> dict[str, StepStatus]:
    result: dict[str, StepStatus] = {}
    for step in steps:
        name = str(step.get("name") or "").strip()
        if not name:
            continue
        result[name] = StepStatus(
            name=name,
            status=normalize(step.get("status")),
            conclusion=normalize(step.get("conclusion")),
        )
    return result


def parse_security_summary(path: Path) -> tuple[dict[str, int] | None, dict[str, int] | None, list[str]]:
    notes: list[str] = []
    try:
        content = path.read_text(encoding="utf-8")
    except OSError as exc:
        return None, None, [f"failed to read summary file: {exc}"]

    bandit_match = BANDIT_RE.search(content)
    pip_match = PIP_AUDIT_RE.search(content)

    bandit: dict[str, int] | None = None
    pip_audit: dict[str, int] | None = None

    if bandit_match:
        bandit = {
            "high": int(bandit_match.group("high")),
            "medium": int(bandit_match.group("medium")),
            "low": int(bandit_match.group("low")),
            "total": int(bandit_match.group("total")),
        }
    else:
        notes.append(
            "bandit summary parse failed: expected "
            "'- Bandit: high=<n>, medium=<n>, low=<n>, total=<n>'"
        )

    if pip_match:
        pip_audit = {
            "vulnerable_packages": int(pip_match.group("pkg")),
            "vulnerabilities": int(pip_match.group("vuln")),
        }
    else:
        notes.append(
            "pip-audit summary parse failed: expected "
            "'- Pip-audit: vulnerable_packages=<n>, vulnerabilities=<n>'"
        )

    return bandit, pip_audit, notes


def find_summary_file(root: Path) -> Path | None:
    for candidate in root.rglob(".ci-security-summary.md"):
        if candidate.is_file():
            return candidate
    return None


def collect_run_evidence(
    run: dict[str, Any],
    repo: str,
    artifact_prefix: str,
) -> RunEvidence:
    db_id = run.get("databaseId")
    if db_id is None:
        raise GhError(f"missing databaseId in run payload: {run.get('id') or run}")
    run_id = int(db_id)
    run_number = int(run.get("number", 0) or 0)
    notes: list[str] = []

    view = run_gh_json(["run", "view", str(run_id), "--repo", repo, "--json", "jobs"])
    jobs = view.get("jobs") or []
    api_job = next((job for job in jobs if str(job.get("name", "")).strip() == "api"), None)

    bandit_step_state = "missing"
    pip_step_state = "missing"
    api_job_conclusion: str | None = None

    if api_job:
        api_job_conclusion = normalize(api_job.get("conclusion"))
        steps = step_map(api_job.get("steps") or [])

        bandit_step = steps.get("Bandit (phase-1 non-blocking)")
        if bandit_step:
            bandit_step_state = bandit_step.conclusion

        pip_step = steps.get("Pip-audit (phase-1 non-blocking)")
        if pip_step:
            pip_step_state = pip_step.conclusion
    else:
        notes.append("api job missing")

    bandit: dict[str, int] | None = None
    pip_audit: dict[str, int] | None = None
    artifact_found = False

    with tempfile.TemporaryDirectory(prefix=f"security-phase1-{run_id}-") as tmp_dir:
        run_dir = Path(tmp_dir) / str(run_id)
        run_dir.mkdir(parents=True, exist_ok=True)
        artifact_name = f"{artifact_prefix}-{run_id}"
        try:
            run_gh(
                [
                    "run",
                    "download",
                    str(run_id),
                    "--repo",
                    repo,
                    "--name",
                    artifact_name,
                    "--dir",
                    str(run_dir),
                ]
            )
            summary_file = find_summary_file(run_dir)
            if summary_file:
                artifact_found = True
                parsed_bandit, parsed_pip, parse_notes = parse_security_summary(summary_file)
                bandit = parsed_bandit
                pip_audit = parsed_pip
                notes.extend(parse_notes)
            else:
                notes.append("security summary artifact missing")
        except GhError as exc:
            notes.append(f"artifact download failed: {exc}")

    baseline_pass: bool | None
    if bandit is not None and pip_audit is not None:
        baseline_pass = bandit["high"] == 0 and pip_audit["vulnerabilities"] == 0
    else:
        baseline_pass = None

    return RunEvidence(
        run_id=run_id,
        run_number=run_number,
        created_at=str(run.get("createdAt") or ""),
        url=str(run.get("url") or ""),
        status=normalize(run.get("status")),
        conclusion=normalize(run.get("conclusion")),
        branch=str(run.get("headBranch") or ""),
        event=str(run.get("event") or ""),
        api_job_conclusion=api_job_conclusion,
        bandit_step=bandit_step_state,
        pip_audit_step=pip_step_state,
        bandit_high=None if bandit is None else bandit["high"],
        bandit_medium=None if bandit is None else bandit["medium"],
        bandit_low=None if bandit is None else bandit["low"],
        bandit_total=None if bandit is None else bandit["total"],
        pip_audit_vulnerable_packages=None if pip_audit is None else pip_audit["vulnerable_packages"],
        pip_audit_vulnerabilities=None if pip_audit is None else pip_audit["vulnerabilities"],
        baseline_pass=baseline_pass,
        artifact_found=artifact_found,
        notes=notes,
    )


def compute_streak(runs: list[RunEvidence], streak_size: int) -> tuple[int, bool]:
    streak = 0
    for run in runs:
        if run.baseline_pass is True:
            streak += 1
            if streak >= streak_size:
                break
        else:
            # Conservative contract: False or unknown (None) both break the latest-pass streak.
            break
    return streak, streak >= streak_size


def iso_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def build_report(
    runs: list[RunEvidence],
    repo: str,
    workflow: str,
    branch: str,
    window_size: int,
    streak_size: int,
) -> dict[str, Any]:
    measured_runs = [run for run in runs if run.baseline_pass is not None]
    pass_count = sum(1 for run in measured_runs if run.baseline_pass)
    fail_count = sum(1 for run in measured_runs if run.baseline_pass is False)
    unknown_count = len(runs) - len(measured_runs)

    streak_count, streak_ok = compute_streak(runs, streak_size)
    window_ok = len(measured_runs) >= window_size
    promotion_ready = window_ok and streak_ok

    return {
        "generated_at": iso_now(),
        "repo": repo,
        "workflow": workflow,
        "branch": branch,
        "criteria": {
            "window_size": window_size,
            "streak_size": streak_size,
            "window_measured_ok": window_ok,
            "latest_streak_ok": streak_ok,
        },
        "stats": {
            "runs_considered": len(runs),
            "runs_measured": len(measured_runs),
            "baseline_pass": pass_count,
            "baseline_fail": fail_count,
            "baseline_unknown": unknown_count,
            "latest_pass_streak": streak_count,
        },
        "promotion_ready": promotion_ready,
        "runs": [asdict(run) for run in runs],
    }


def format_markdown(report: dict[str, Any]) -> str:
    criteria = report["criteria"]
    stats = report["stats"]
    lines = [
        "## Security Phase-1 Promotion Evidence",
        "",
        f"- Repo: `{report['repo']}`",
        f"- Workflow: `{report['workflow']}`",
        f"- Branch: `{report['branch']}`",
        f"- Generated (UTC): `{report['generated_at']}`",
        f"- Promotion ready: `{'yes' if report['promotion_ready'] else 'no'}`",
        "",
        "### Criteria status",
        "",
        f"- {criteria['window_size']}-run window measured: `{criteria['window_measured_ok']}` ({stats['runs_measured']}/{criteria['window_size']})",
        f"- Latest {criteria['streak_size']} consecutive baseline stability: `{criteria['latest_streak_ok']}` (streak={stats['latest_pass_streak']}/{criteria['streak_size']})",
        "",
        "### Aggregate counts",
        "",
        f"- Baseline pass: `{stats['baseline_pass']}`",
        f"- Baseline fail: `{stats['baseline_fail']}`",
        f"- Baseline unknown: `{stats['baseline_unknown']}`",
        "",
        "### Recent runs",
        "",
        "| Run | Created | Conclusion | Bandit high | Pip-audit vulns | Baseline | Notes |",
        "| --- | --- | --- | ---: | ---: | --- | --- |",
    ]

    for run in report["runs"]:
        run_link = f"[{run['run_number']}]({run['url']})" if run["url"] else str(run["run_number"])
        notes = "; ".join(run["notes"]) if run["notes"] else ""
        lines.append(
            "| "
            + " | ".join(
                [
                    run_link,
                    run["created_at"] or "",
                    run["conclusion"],
                    "" if run["bandit_high"] is None else str(run["bandit_high"]),
                    "" if run["pip_audit_vulnerabilities"] is None else str(run["pip_audit_vulnerabilities"]),
                    "pass" if run["baseline_pass"] is True else "fail" if run["baseline_pass"] is False else "unknown",
                    notes.replace("|", "\\|"),
                ]
            )
            + " |"
        )

    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate security phase-1 promotion evidence from GitHub runs")
    parser.add_argument("--repo", default=None, help="Repository in OWNER/REPO format (defaults to current repo)")
    parser.add_argument("--workflow", default="ci.yml", help="Workflow name or file (default: ci.yml)")
    parser.add_argument("--branch", default="main", help="Target branch (default: main)")
    parser.add_argument("--event", default="push", help="Event filter for runs (default: push)")
    parser.add_argument("--window-size", type=int, default=10, help="Required measured run window")
    parser.add_argument("--streak-size", type=int, default=5, help="Required latest pass streak")
    parser.add_argument("--run-limit", type=int, default=20, help="Max runs to inspect (must be >= window-size)")
    parser.add_argument(
        "--artifact-prefix",
        default="api-security-phase1-reports",
        help="Artifact name prefix used in CI upload",
    )
    parser.add_argument("--markdown-out", default=None, help="Optional path to write markdown report")
    parser.add_argument("--json-out", default=None, help="Optional path to write JSON report")
    parser.add_argument(
        "--format",
        choices=["markdown", "json", "both"],
        default="markdown",
        help="Primary stdout format",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero if promotion criteria are not met",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.window_size <= 0 or args.streak_size <= 0:
        print("window-size and streak-size must be > 0", file=sys.stderr)
        return 3
    if args.run_limit < args.window_size:
        print("run-limit must be >= window-size", file=sys.stderr)
        return 3

    try:
        repo = args.repo
        if not repo:
            repo_payload = run_gh_json(["repo", "view", "--json", "nameWithOwner"])
            repo = str(repo_payload["nameWithOwner"])

        runs_payload = run_gh_json(
            [
                "run",
                "list",
                "--repo",
                repo,
                "--workflow",
                args.workflow,
                "--branch",
                args.branch,
                "--event",
                args.event,
                "--limit",
                str(args.run_limit),
                "--json",
                "databaseId,number,createdAt,url,status,conclusion,headBranch,event",
            ]
        )
        if not isinstance(runs_payload, list):
            raise GhError("unexpected response for gh run list")

        # Keep this serial for now to limit scope; optimize N+1 gh calls separately if needed.
        run_evidence = [
            collect_run_evidence(run, repo=repo, artifact_prefix=args.artifact_prefix)
            for run in runs_payload
        ]
        report = build_report(
            runs=run_evidence,
            repo=repo,
            workflow=args.workflow,
            branch=args.branch,
            window_size=args.window_size,
            streak_size=args.streak_size,
        )

        markdown_text = format_markdown(report)
        json_text = json.dumps(report, indent=2, sort_keys=True) + "\n"

        if args.markdown_out:
            md_path = Path(args.markdown_out)
            md_path.parent.mkdir(parents=True, exist_ok=True)
            md_path.write_text(markdown_text, encoding="utf-8")

        if args.json_out:
            json_path = Path(args.json_out)
            json_path.parent.mkdir(parents=True, exist_ok=True)
            json_path.write_text(json_text, encoding="utf-8")

        if args.format == "markdown":
            print(markdown_text, end="")
        elif args.format == "json":
            print(json_text, end="")
        else:
            print(markdown_text, end="")
            print("\n```json")
            print(json_text, end="")
            print("```")

        if args.strict and not report["promotion_ready"]:
            return 2
        return 0

    except Exception as exc:  # noqa: BLE001
        print(str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
