#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Guard CI_SELF_HOSTED_TRUSTED_ROUTING with local/GitHub runner health.

Usage:
  scripts/runner/guard_trusted_routing.sh \
    --owner grace-bidos \
    --repo chem-model-edit \
    --labels "self-hosted,linux,x64,chem-trusted-pr"

Options:
  --owner <owner>                Required.
  --repo <repo>                  Required.
  --labels <csv>                 Optional. Default: self-hosted,linux,x64,chem-trusted-pr
  --var-name <name>              Optional. Default: CI_SELF_HOSTED_TRUSTED_ROUTING
  --no-toggle                    Do not write the variable; health check only.
  --strict-gh                    Pass --strict-gh to health checker.
  --allow-zero-capacity          Pass --allow-zero-capacity to health checker.
  --dry-run                      Print actions only.
  -h, --help                     Show help.

Behavior:
  - Runs check_local_runner_health.sh.
  - If health is degraded and toggle is enabled, sets variable to false.
  - If health is healthy, leaves variable unchanged.
USAGE
}

die() {
  echo "$*" >&2
  exit 1
}

owner=""
repo=""
labels="self-hosted,linux,x64,chem-trusted-pr"
var_name="CI_SELF_HOSTED_TRUSTED_ROUTING"
no_toggle=0
strict_gh=0
allow_zero_capacity=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --owner) [[ $# -ge 2 ]] || die "--owner requires a value"; owner="${2:-}"; shift 2 ;;
    --repo) [[ $# -ge 2 ]] || die "--repo requires a value"; repo="${2:-}"; shift 2 ;;
    --labels) [[ $# -ge 2 ]] || die "--labels requires a value"; labels="${2:-}"; shift 2 ;;
    --var-name) [[ $# -ge 2 ]] || die "--var-name requires a value"; var_name="${2:-}"; shift 2 ;;
    --no-toggle) no_toggle=1; shift ;;
    --strict-gh) strict_gh=1; shift ;;
    --allow-zero-capacity) allow_zero_capacity=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "$owner" ]] || die "--owner is required"
[[ -n "$repo" ]] || die "--repo is required"

for tool in gh; do
  command -v "$tool" >/dev/null 2>&1 || die "missing command: $tool"
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
health_script="${script_dir}/check_local_runner_health.sh"
[[ -x "$health_script" ]] || die "missing executable: $health_script"

health_args=(
  --owner "$owner"
  --repo "$repo"
  --labels "$labels"
)
if [[ "$strict_gh" -eq 1 ]]; then
  health_args+=(--strict-gh)
fi
if [[ "$allow_zero_capacity" -eq 1 ]]; then
  health_args+=(--allow-zero-capacity)
fi

echo "Running health check for ${owner}/${repo}"
if "$health_script" "${health_args[@]}"; then
  echo "Health is OK. No routing change applied."
  exit 0
fi

echo "Health is DEGRADED."
if [[ "$no_toggle" -eq 1 ]]; then
  echo "--no-toggle set; not changing ${var_name}."
  exit 1
fi

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] gh variable set ${var_name} --repo ${owner}/${repo} --body false"
  exit 1
fi

gh variable set "$var_name" --repo "${owner}/${repo}" --body false

echo "Set ${var_name}=false for ${owner}/${repo}."
exit 1
