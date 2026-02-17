#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
One-command setup for always-on ephemeral runner pool supervisor.

Usage:
  scripts/runner/setup_pool_supervisor_one_command.sh [options]

Options:
  --repo <owner/repo>               Default: grace-bidos/chem-model-edit
  --min <n>                         Default: 1
  --baseline <n>                    Deprecated alias for --min
  --max <n>                         Default: 4
  --target <n>                      Default: max
  --interval <seconds>              Default: 15
  --lock-file <path>                Default: /run/lock/chem-model-edit-runner-pool.lock
  --gh-token-file <path>            Default: /etc/chem-model-edit/gh_runner_token
  --token-source <gh|env|file|app>  Default: gh
  --app-id <id>                     Required for app mode.
  --app-installation-id <id>        Required for app mode.
  --app-private-key-file <path>     Required for app mode.
  --disable-timer                   Disable old reconcile timer (default: on)
  --dry-run                         Print actions only.
  -h, --help                        Show help.

Token source behavior:
  gh   : use `gh auth token`
  env  : use existing GH_TOKEN env var
  file : do not write token file; only use existing --gh-token-file
  app  : request short-lived installation token via GitHub App credentials
EOF
}

repo="grace-bidos/chem-model-edit"
min_pool="1"
max_parallel="4"
target=""
interval_seconds="15"
lock_file="/run/lock/chem-model-edit-runner-pool.lock"
gh_token_file="/etc/chem-model-edit/gh_runner_token"
token_source="gh"
app_id=""
app_installation_id=""
app_private_key_file=""
disable_timer=1
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
    --gh-token-file) gh_token_file="${2:-}"; shift 2 ;;
    --token-source) token_source="${2:-}"; shift 2 ;;
    --app-id) app_id="${2:-}"; shift 2 ;;
    --app-installation-id) app_installation_id="${2:-}"; shift 2 ;;
    --app-private-key-file) app_private_key_file="${2:-}"; shift 2 ;;
    --disable-timer) disable_timer=1; shift ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$target" ]]; then
  target="$max_parallel"
fi

if [[ "$token_source" != "gh" && "$token_source" != "env" && "$token_source" != "file" && "$token_source" != "app" ]]; then
  echo "--token-source must be one of: gh, env, file, app" >&2
  exit 1
fi

for tool in bash sudo; do
  command -v "$tool" >/dev/null 2>&1 || { echo "missing command: $tool" >&2; exit 1; }
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
install_script="${script_dir}/install_pool_supervisor_service.sh"
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
elif [[ "$token_source" == "app" ]]; then
  if [[ -z "$app_id" || -z "$app_installation_id" || -z "$app_private_key_file" ]]; then
    echo "app mode requires --app-id, --app-installation-id, and --app-private-key-file" >&2
    exit 1
  fi
  app_token_script="${script_dir}/request_github_app_token.sh"
  if [[ ! -x "$app_token_script" ]]; then
    echo "missing executable: $app_token_script" >&2
    exit 1
  fi
  token_value="$("$app_token_script" \
    --app-id "$app_id" \
    --installation-id "$app_installation_id" \
    --private-key-file "$app_private_key_file")"
  if [[ -z "$token_value" ]]; then
    echo "Failed to obtain GitHub App installation token." >&2
    exit 1
  fi
fi

if [[ "$dry_run" -eq 1 ]]; then
  echo "[dry-run] repo=${repo} min=${min_pool} max=${max_parallel} target=${target} interval=${interval_seconds}"
  echo "[dry-run] lock_file=${lock_file} token_source=${token_source} gh_token_file=${gh_token_file}"
  if [[ "$token_source" == "gh" || "$token_source" == "env" || "$token_source" == "app" ]]; then
    echo "[dry-run] token_length=$(printf '%s' "$token_value" | wc -c | tr -d ' ')"
    echo "[dry-run] sudo install -d -m 0700 $(dirname "$gh_token_file")"
    echo "[dry-run] sudo sh -c \"umask 077; printf '%s\\n' '***' > '$gh_token_file'\""
  else
    echo "[dry-run] token file must already exist: $gh_token_file"
  fi
  "$install_script" \
    --repo "$repo" \
    --min "$min_pool" \
    --max "$max_parallel" \
    --target "$target" \
    --interval "$interval_seconds" \
    --lock-file "$lock_file" \
    --gh-token-file "$gh_token_file" \
    --disable-timer \
    --dry-run
  exit 0
fi

if ! sudo -n true >/dev/null 2>&1; then
  echo "sudo auth is required. Run: sudo -v" >&2
  exit 1
fi

if [[ "$token_source" == "gh" || "$token_source" == "env" || "$token_source" == "app" ]]; then
  sudo install -d -m 0700 "$(dirname "$gh_token_file")"
  sudo sh -c "umask 077; printf '%s\n' '$token_value' > '$gh_token_file'"
fi

install_args=(
  --repo "$repo"
  --min "$min_pool"
  --max "$max_parallel"
  --target "$target"
  --interval "$interval_seconds"
  --lock-file "$lock_file"
  --gh-token-file "$gh_token_file"
)
if [[ "$disable_timer" -eq 1 ]]; then
  install_args+=(--disable-timer)
fi

"$install_script" "${install_args[@]}"

sudo systemctl status chem-runner-pool-supervisor.service --no-pager
