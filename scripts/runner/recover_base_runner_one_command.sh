#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  scripts/runner/recover_base_runner_one_command.sh [options]

Options:
  --owner <org-or-user>        Default: grace-bidos
  --repo <repo-name>           Default: chem-model-edit
  --labels <csv>               Default: self-hosted,linux,x64,chem-trusted-pr
  --group <runner-group>       Default: Default
  --name <runner-name>         Default: home-self-host
  --service-user <user>        Default: runner-user
  --runner-home <path>         Default: /opt/actions-runner/actions-runner
  --dry-run                    Show actions without mutating runner state
  -h, --help                   Show this help

What this wrapper does:
  1) Ensures required tools exist (`gh`, `sudo`).
  2) Fetches GH token from `gh auth token` if `GH_TOKEN` is unset.
  3) Runs `recover_base_runner.sh` with required envs.
     (`recover_base_runner.sh` itself uses sudo only for service operations)
USAGE
}

owner="grace-bidos"
repo="chem-model-edit"
labels="self-hosted,linux,x64,chem-trusted-pr"
runner_group="Default"
runner_name="home-self-host"
service_user="runner-user"
runner_home="/opt/actions-runner/actions-runner"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --owner) owner="${2:-}"; shift 2 ;;
    --repo) repo="${2:-}"; shift 2 ;;
    --labels) labels="${2:-}"; shift 2 ;;
    --group) runner_group="${2:-}"; shift 2 ;;
    --name) runner_name="${2:-}"; shift 2 ;;
    --service-user) service_user="${2:-}"; shift 2 ;;
    --runner-home) runner_home="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

for tool in gh sudo; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing required tool: $tool" >&2
    exit 2
  fi
done

if [[ -z "${GH_TOKEN:-}" ]]; then
  if ! gh auth status -h github.com >/dev/null 2>&1; then
    echo "gh is not authenticated for github.com. Run: gh auth login" >&2
    exit 1
  fi
  GH_TOKEN="$(gh auth token)"
fi

if [[ -z "${GH_TOKEN}" ]]; then
  echo "Failed to obtain GH_TOKEN." >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
recover_script="${script_dir}/recover_base_runner.sh"

if [[ ! -x "${recover_script}" ]]; then
  echo "Missing executable: ${recover_script}" >&2
  exit 1
fi

export GH_TOKEN
export RUNNER_OWNER="${owner}"
export RUNNER_REPO="${repo}"
export RUNNER_LABELS="${labels}"
export RUNNER_GROUP="${runner_group}"
export RUNNER_NAME="${runner_name}"
export RUNNER_SERVICE_USER="${service_user}"

cmd=(
  "${recover_script}"
  "--runner-home" "${runner_home}"
)

if [[ "${dry_run}" -eq 1 ]]; then
  cmd+=("--dry-run")
fi

echo "Running base runner recovery with:"
echo "  owner=${owner}"
echo "  repo=${repo}"
echo "  labels=${labels}"
echo "  group=${runner_group}"
echo "  name=${runner_name}"
echo "  service_user=${service_user}"
echo "  runner_home=${runner_home}"
if [[ "${dry_run}" -eq 1 ]]; then
  echo "  mode=dry-run"
fi

"${cmd[@]}"
