#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
One-command setup for ephemeral runner auto-reconcile timer.

Usage:
  scripts/runner/setup_pool_reconcile_one_command.sh [options]

Options:
  --repo <owner/repo>               Default: grace-bidos/chem-model-edit
  --baseline <n>                    Default: 1
  --max <n>                         Default: 4
  --target <n>                      Default: max
  --interval <duration>             Default: 2min
  --gh-token-file <path>            Default: /etc/chem-model-edit/gh_runner_token
  --token-source <gh|env|file>      Default: gh
  --dry-run                         Print actions only
  -h, --help                        Show help

Token source behavior:
  gh   : use `gh auth token`
  env  : use existing GH_TOKEN env var
  file : do not write token file; only use existing --gh-token-file
EOF
}

repo="grace-bidos/chem-model-edit"
baseline="1"
max_parallel="4"
target=""
interval="2min"
gh_token_file="/etc/chem-model-edit/gh_runner_token"
token_source="gh"
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo) repo="${2:-}"; shift 2 ;;
    --baseline) baseline="${2:-}"; shift 2 ;;
    --max) max_parallel="${2:-}"; shift 2 ;;
    --target) target="${2:-}"; shift 2 ;;
    --interval) interval="${2:-}"; shift 2 ;;
    --gh-token-file) gh_token_file="${2:-}"; shift 2 ;;
    --token-source) token_source="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$target" ]]; then
  target="$max_parallel"
fi

if [[ "$token_source" != "gh" && "$token_source" != "env" && "$token_source" != "file" ]]; then
  echo "--token-source must be one of: gh, env, file" >&2
  exit 1
fi

for tool in bash sudo; do
  command -v "$tool" >/dev/null 2>&1 || { echo "missing command: $tool" >&2; exit 1; }
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
install_script="${script_dir}/install_pool_reconcile_timer.sh"
if [[ ! -x "$install_script" ]]; then
  echo "missing executable: $install_script" >&2
  exit 1
fi

token_value=""
if [[ "$token_source" == "gh" ]]; then
  command -v gh >/dev/null 2>&1 || { echo "missing command: gh" >&2; exit 1; }
  if ! gh auth status -h github.com >/dev/null 2>&1; then
    echo "gh is not authenticated. Run: gh auth login" >&2
    exit 1
  fi
  token_value="$(gh auth token)"
elif [[ "$token_source" == "env" ]]; then
  token_value="${GH_TOKEN:-}"
  if [[ -z "$token_value" ]]; then
    echo "GH_TOKEN is empty with --token-source env" >&2
    exit 1
  fi
fi

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] repo=${repo} baseline=${baseline} max=${max_parallel} target=${target} interval=${interval}"
  echo "[dry-run] token_source=${token_source} gh_token_file=${gh_token_file}"
  if [[ "$token_source" == "gh" || "$token_source" == "env" ]]; then
    echo "[dry-run] token_length=$(printf '%s' "$token_value" | wc -c | tr -d ' ')"
    echo "[dry-run] sudo install -d -m 0700 $(dirname "$gh_token_file")"
    echo "[dry-run] sudo sh -c \"umask 077; printf '%s\\n' '***' > '$gh_token_file'\""
  else
    echo "[dry-run] token file must already exist: $gh_token_file"
  fi
  echo "[dry-run] sudo systemctl disable --now chem-runner-pool-reconcile.timer || true"
  echo "[dry-run] sudo systemctl stop chem-runner-pool-reconcile.service || true"
  "$install_script" \
    --repo "$repo" \
    --baseline "$baseline" \
    --max "$max_parallel" \
    --target "$target" \
    --interval "$interval" \
    --gh-token-file "$gh_token_file" \
    --dry-run
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

if [[ "$token_source" == "gh" || "$token_source" == "env" ]]; then
  sudo install -d -m 0700 "$(dirname "$gh_token_file")"
  sudo sh -c "umask 077; printf '%s\n' '$token_value' > '$gh_token_file'"
fi

sudo systemctl disable --now chem-runner-pool-reconcile.timer >/dev/null 2>&1 || true
sudo systemctl stop chem-runner-pool-reconcile.service >/dev/null 2>&1 || true

"$install_script" \
  --repo "$repo" \
  --baseline "$baseline" \
  --max "$max_parallel" \
  --target "$target" \
  --interval "$interval" \
  --gh-token-file "$gh_token_file"

sudo systemctl status chem-runner-pool-reconcile.timer --no-pager
sudo systemctl status chem-runner-pool-reconcile.service --no-pager || true
