from __future__ import annotations

import json
from pathlib import Path
import os
import stat
import subprocess

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = REPO_ROOT / "scripts" / "runner" / "scale_local_runner_pool.sh"


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _prepare_stub_bin(tmp_path: Path) -> Path:
    bindir = tmp_path / "bin"
    bindir.mkdir(parents=True, exist_ok=True)

    _write_executable(
        bindir / "gh",
        """#!/usr/bin/env bash
set -euo pipefail

if [[ "${1:-}" != "api" ]]; then
  echo "unsupported gh stub command" >&2
  exit 1
fi
shift

method="GET"
if [[ "${1:-}" == "-X" ]]; then
  method="${2:-GET}"
  shift 2
fi

endpoint="${1:-}"
shift || true

if [[ "$method" == "DELETE" ]]; then
  if [[ -n "${GH_STUB_LOG:-}" ]]; then
    printf 'DELETE %s\n' "$endpoint" >>"$GH_STUB_LOG"
  fi
  echo '{}'
  exit 0
fi

if [[ "$endpoint" == repos/*/actions/runners/registration-token ]]; then
  echo '{"token":"stub-registration-token"}'
  exit 0
fi

if [[ "$endpoint" == repos/*/actions/runners/remove-token ]]; then
  echo '{"token":"stub-remove-token"}'
  exit 0
fi

if [[ "$endpoint" == repos/*/actions/runners ]]; then
  payload="${GH_STUB_RUNNERS_JSON:-{\"runners\":[]}}"
  if [[ "${1:-}" == "--jq" ]]; then
    jq -r "${2:-.}" <<<"$payload"
  else
    echo "$payload"
  fi
  exit 0
fi

echo '{}'
""",
    )

    for tool in ("systemctl", "rsync", "sudo", "runuser"):
        _write_executable(
            bindir / tool,
            "#!/usr/bin/env bash\nset -euo pipefail\nexit 0\n",
        )

    return bindir


def _run_scaler(args: list[str], env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(SCRIPT_PATH), *args],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )


def test_scale_local_runner_pool_emits_machine_readable_inventory_and_status(tmp_path: Path) -> None:
    bindir = _prepare_stub_bin(tmp_path)
    base_dir = tmp_path / "runner-base"
    runner1 = base_dir / "runner"
    runner1.mkdir(parents=True)
    _write_executable(runner1 / "config.sh", "#!/usr/bin/env bash\nexit 0\n")
    (runner1 / ".runner").write_text("configured", encoding="utf-8")

    inventory_path = tmp_path / "inventory.json"
    status_path = tmp_path / "status.json"

    gh_payload = {
        "runners": [
            {"id": 104, "name": "runner-4", "busy": False, "status": "online"},
            {"id": 103, "name": "runner-3", "busy": False, "status": "online"},
            {"id": 102, "name": "runner-2", "busy": False, "status": "online"},
            {"id": 101, "name": "runner", "busy": False, "status": "online"},
        ]
    }

    env = os.environ.copy()
    env["PATH"] = f"{bindir}:{env['PATH']}"
    env["GH_STUB_RUNNERS_JSON"] = json.dumps(gh_payload)

    result = _run_scaler(
        [
            "--repo",
            "grace-bidos/chem-model-edit",
            "--base-dir",
            str(base_dir),
            "--base-name",
            "runner",
            "--min",
            "2",
            "--max",
            "4",
            "--target",
            "1",
            "--dry-run",
            "--inventory-out",
            str(inventory_path),
            "--status-out",
            str(status_path),
        ],
        env,
    )

    assert result.returncode == 0, result.stderr
    assert "min=2 max=4 requested_target=1 effective_target=2" in result.stdout
    assert result.stdout.index("remove: runner-4") < result.stdout.index("remove: runner-3")

    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))

    assert inventory["control"] == {
        "min": 2,
        "max": 4,
        "target_requested": 1,
        "target_effective": 2,
    }
    assert status["counts"]["warnings"] == 0
    assert status["counts"]["operations"] == 4

    operation_actions = [entry["action"] for entry in status["operations"]]
    assert operation_actions == ["ensure_started", "create", "remove", "remove"]


def test_scale_local_runner_pool_rejects_conflicting_min_and_baseline(tmp_path: Path) -> None:
    bindir = _prepare_stub_bin(tmp_path)
    env = os.environ.copy()
    env["PATH"] = f"{bindir}:{env['PATH']}"
    env["GH_STUB_RUNNERS_JSON"] = '{"runners":[]}'

    result = _run_scaler(
        [
            "--repo",
            "grace-bidos/chem-model-edit",
            "--min",
            "2",
            "--baseline",
            "3",
            "--dry-run",
        ],
        env,
    )

    assert result.returncode == 1
    assert "--min and --baseline conflict" in result.stderr
