#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Install systemd timer/service to keep ephemeral runner pool replenished.

Usage:
  scripts/runner/install_pool_reconcile_timer.sh \
    --repo grace-bidos/chem-model-edit \
    --min 1 \
    --max 4 \
    --target 4 \
    --interval 2min \
    --gh-token-file /etc/chem-model-edit/gh_runner_token

Options:
  --repo <owner/repo>               Required.
  --min <n>                         Optional. Default: 1
  --baseline <n>                    Deprecated alias for --min.
  --max <n>                         Optional. Default: 4
  --target <n>                      Optional. Default: max
  --lock-file <path>                Optional. Default: /run/lock/chem-model-edit-runner-pool.lock
  --interval <systemd duration>     Optional. Default: 2min
  --service-name <name>             Optional. Default: chem-runner-pool-reconcile
  --gh-token-file <path>            Required. Root-readable file with GH token.
  --repo-root <path>                Optional. Default: current git root
  --dry-run                         Print files/commands only.
  -h, --help                        Show help.
EOF
}

repo=""
min_pool="1"
max_parallel="4"
target=""
lock_file="/run/lock/chem-model-edit-runner-pool.lock"
interval="2min"
service_name="chem-runner-pool-reconcile"
gh_token_file=""
repo_root=""
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) repo="${2:-}"; shift 2 ;;
    --min) min_pool="${2:-}"; shift 2 ;;
    --baseline) min_pool="${2:-}"; shift 2 ;;
    --max) max_parallel="${2:-}"; shift 2 ;;
    --target) target="${2:-}"; shift 2 ;;
    --lock-file) lock_file="${2:-}"; shift 2 ;;
    --interval) interval="${2:-}"; shift 2 ;;
    --service-name) service_name="${2:-}"; shift 2 ;;
    --gh-token-file) gh_token_file="${2:-}"; shift 2 ;;
    --repo-root) repo_root="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$repo" || -z "$gh_token_file" ]]; then
  echo "--repo and --gh-token-file are required." >&2
  usage
  exit 1
fi

if [[ -z "$repo_root" ]]; then
  repo_root="$(git rev-parse --show-toplevel)"
fi
if [[ -z "$target" ]]; then
  target="$max_parallel"
fi

reconcile_script="${repo_root}/scripts/runner/reconcile_ephemeral_pool.sh"
if [[ ! -x "$reconcile_script" ]]; then
  echo "missing executable: $reconcile_script" >&2
  exit 1
fi

service_file="/etc/systemd/system/${service_name}.service"
timer_file="/etc/systemd/system/${service_name}.timer"

service_content="[Unit]
Description=Reconcile local GitHub ephemeral runner pool
After=network-online.target

[Service]
Type=oneshot
ExecStart=/usr/bin/env bash -lc 'set -euo pipefail; token_file=${gh_token_file}; export GH_TOKEN=\"\$(cat \"\$token_file\")\"; ${reconcile_script} --repo ${repo} --min ${min_pool} --max ${max_parallel} --target ${target} --lock-file ${lock_file}'
"

timer_content="[Unit]
Description=Periodic reconcile of local GitHub ephemeral runner pool

[Timer]
OnBootSec=30s
OnUnitActiveSec=${interval}
Persistent=true

[Install]
WantedBy=timers.target
"

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] write ${service_file}:"
  printf '%s\n' "$service_content"
  echo "[dry-run] write ${timer_file}:"
  printf '%s\n' "$timer_content"
  echo "[dry-run] sudo install token file permissions check"
  echo "[dry-run] sudo systemctl daemon-reload"
  echo "[dry-run] sudo systemctl enable --now ${service_name}.timer"
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

tmpd="$(mktemp -d)"
trap 'rm -rf "$tmpd"' EXIT

printf '%s\n' "$service_content" > "${tmpd}/service"
printf '%s\n' "$timer_content" > "${tmpd}/timer"

sudo install -m 0644 "${tmpd}/service" "$service_file"
sudo install -m 0644 "${tmpd}/timer" "$timer_file"

sudo systemctl daemon-reload
sudo systemctl enable --now "${service_name}.timer"
sudo systemctl start "${service_name}.service"
sudo systemctl status "${service_name}.timer" --no-pager
