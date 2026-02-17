#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Install systemd supervisor service to keep ephemeral runner pool at target size.

Usage:
  scripts/runner/install_pool_supervisor_service.sh \
    --repo grace-bidos/chem-model-edit \
    --min 1 \
    --max 4 \
    --target 4 \
    --interval 15 \
    --gh-token-file /etc/chem-model-edit/gh_runner_token

Options:
  --repo <owner/repo>               Required.
  --min <n>                         Optional. Default: 1
  --baseline <n>                    Deprecated alias for --min.
  --max <n>                         Optional. Default: 4
  --target <n>                      Optional. Default: max
  --interval <seconds>              Optional. Default: 15
  --lock-file <path>                Optional. Default: /run/lock/chem-model-edit-runner-pool.lock
  --service-name <name>             Optional. Default: chem-runner-pool-supervisor
  --gh-token-file <path>            Required. Root-readable file with GH token.
  --repo-root <path>                Optional. Default: current git root
  --disable-timer                   Disable old reconcile timer after install.
  --dry-run                         Print files/commands only.
  -h, --help                        Show help.
EOF
}

repo=""
min_pool="1"
max_parallel="4"
target=""
interval_seconds="15"
lock_file="/run/lock/chem-model-edit-runner-pool.lock"
service_name="chem-runner-pool-supervisor"
gh_token_file=""
repo_root=""
disable_timer=0
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) repo="${2:-}"; shift 2 ;;
    --min) min_pool="${2:-}"; shift 2 ;;
    --baseline) min_pool="${2:-}"; shift 2 ;;
    --max) max_parallel="${2:-}"; shift 2 ;;
    --target) target="${2:-}"; shift 2 ;;
    --interval) interval_seconds="${2:-}"; shift 2 ;;
    --lock-file) lock_file="${2:-}"; shift 2 ;;
    --service-name) service_name="${2:-}"; shift 2 ;;
    --gh-token-file) gh_token_file="${2:-}"; shift 2 ;;
    --repo-root) repo_root="${2:-}"; shift 2 ;;
    --disable-timer) disable_timer=1; shift ;;
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
if [[ -z "$target" ]]; then
  target="$max_parallel"
fi
if [[ -z "$repo_root" ]]; then
  repo_root="$(git rev-parse --show-toplevel)"
fi

supervisor_script="${repo_root}/scripts/runner/supervise_ephemeral_pool.sh"
if [[ ! -x "$supervisor_script" ]]; then
  echo "missing executable: $supervisor_script" >&2
  exit 1
fi

service_file="/etc/systemd/system/${service_name}.service"

service_content="[Unit]
Description=Supervisor for local GitHub ephemeral runner pool
After=network-online.target
Wants=network-online.target

[Service]
Type=simple
ExecStart=/usr/bin/env bash -lc 'set -euo pipefail; token_file=${gh_token_file}; export GH_TOKEN=\"\$(cat \"\$token_file\")\"; ${supervisor_script} --repo ${repo} --min ${min_pool} --max ${max_parallel} --target ${target} --interval ${interval_seconds} --lock-file ${lock_file}'
Restart=always
RestartSec=5s
KillMode=process

[Install]
WantedBy=multi-user.target
"

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] write ${service_file}:"
  printf '%s\n' "$service_content"
  echo "[dry-run] sudo systemctl daemon-reload"
  echo "[dry-run] sudo systemctl enable --now ${service_name}.service"
  if [[ "$disable_timer" -eq 1 ]]; then
    echo "[dry-run] sudo systemctl disable --now chem-runner-pool-reconcile.timer || true"
    echo "[dry-run] sudo systemctl stop chem-runner-pool-reconcile.service || true"
  fi
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

tmpd="$(mktemp -d)"
trap 'rm -rf "$tmpd"' EXIT
printf '%s\n' "$service_content" > "${tmpd}/service"

sudo install -m 0644 "${tmpd}/service" "$service_file"
sudo systemctl daemon-reload
sudo systemctl enable --now "${service_name}.service"

if [[ "$disable_timer" -eq 1 ]]; then
  sudo systemctl disable --now chem-runner-pool-reconcile.timer >/dev/null 2>&1 || true
  sudo systemctl stop chem-runner-pool-reconcile.service >/dev/null 2>&1 || true
fi

sudo systemctl status "${service_name}.service" --no-pager
