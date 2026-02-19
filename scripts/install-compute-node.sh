#!/usr/bin/env bash
set -euo pipefail

API_BASE="http://localhost:8000"
JOIN_TOKEN=""
NODE_NAME="$(hostname -s 2>/dev/null || hostname || echo compute-node)"
QUEUE_NAME=""

usage() {
  cat <<'USAGE'
Usage: install-compute-node.sh [options]

Register this machine as a runtime compute node using a short-lived join token.

Options:
  --api-base <url>       Runtime API base URL (default: http://localhost:8000)
  --join-token <token>   Required join token from /api/runtime/nodes/join-token
  --name <value>         Logical node name (default: hostname)
  --queue-name <value>   Optional queue hint for metadata
  -h, --help             Show this help
USAGE
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

have_command() {
  command -v "$1" >/dev/null 2>&1
}

normalize_api_base() {
  local value="$1"
  value="${value%/}"
  if [[ "$value" == */api ]]; then
    value="${value%/api}"
  fi
  printf '%s\n' "$value"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --api-base)
      shift
      [[ $# -gt 0 ]] || die "--api-base requires a value"
      API_BASE="$1"
      ;;
    --join-token)
      shift
      [[ $# -gt 0 ]] || die "--join-token requires a value"
      JOIN_TOKEN="$1"
      ;;
    --name)
      shift
      [[ $# -gt 0 ]] || die "--name requires a value"
      NODE_NAME="$1"
      ;;
    --queue-name)
      shift
      [[ $# -gt 0 ]] || die "--queue-name requires a value"
      QUEUE_NAME="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
  shift
done

[[ -n "$JOIN_TOKEN" ]] || die "--join-token is required"
have_command curl || die "curl is required"
have_command jq || die "jq is required"

API_BASE="$(normalize_api_base "$API_BASE")"
REGISTER_URL="$API_BASE/api/runtime/nodes/register"

sinfo_ok=false
sbatch_ok=false
squeue_ok=false
sacct_ok=false
slurmd_active=false
munge_active=false
slurm_state="unknown"

if have_command sinfo; then
  if sinfo --version >/dev/null 2>&1; then
    sinfo_ok=true
  fi
fi

if have_command sbatch; then
  if sbatch --version >/dev/null 2>&1; then
    sbatch_ok=true
  fi
fi

if have_command squeue; then
  if squeue --version >/dev/null 2>&1; then
    squeue_ok=true
  fi
fi

if have_command sacct; then
  if sacct --version >/dev/null 2>&1; then
    sacct_ok=true
  fi
fi

if have_command systemctl; then
  if systemctl is-active --quiet slurmd; then
    slurmd_active=true
  fi
  if systemctl is-active --quiet munge; then
    munge_active=true
  fi
fi

if [[ "$sinfo_ok" == "true" ]]; then
  slurm_state="$(sinfo -h -o %T 2>/dev/null | head -n1 | tr -d '[:space:]')"
  if [[ -z "$slurm_state" ]]; then
    slurm_state="unknown"
  fi
fi

payload="$(jq -n \
  --arg token "$JOIN_TOKEN" \
  --arg name "$NODE_NAME" \
  --arg host "$(hostname -f 2>/dev/null || hostname || echo unknown-host)" \
  --arg queue_hint "$QUEUE_NAME" \
  --arg slurm_state "$slurm_state" \
  --arg os "$(uname -s)" \
  --arg kernel "$(uname -r)" \
  --argjson sinfo "$sinfo_ok" \
  --argjson sbatch "$sbatch_ok" \
  --argjson squeue "$squeue_ok" \
  --argjson sacct "$sacct_ok" \
  --argjson slurmd_active "$slurmd_active" \
  --argjson munge_active "$munge_active" \
  '{
    token: $token,
    name: $name,
    meta: {
      host: $host,
      queue_hint: $queue_hint,
      slurm_state: $slurm_state,
      os: $os,
      kernel: $kernel,
      health_checks: {
        sinfo: $sinfo,
        sbatch: $sbatch,
        squeue: $squeue,
        sacct: $sacct,
        slurmd_active: $slurmd_active,
        munge_active: $munge_active
      }
    }
  }')"

response="$(
  curl -fsS \
    -X POST \
    -H 'content-type: application/json' \
    "$REGISTER_URL" \
    -d "$payload"
)"

printf 'Compute node registered\n'
printf '  server_id: %s\n' "$(printf '%s' "$response" | jq -r '.server_id')"
printf '  target_id: %s\n' "$(printf '%s' "$response" | jq -r '.target_id')"
printf '  queue_name: %s\n' "$(printf '%s' "$response" | jq -r '.queue_name')"
printf '  registered_at: %s\n' "$(printf '%s' "$response" | jq -r '.registered_at')"
