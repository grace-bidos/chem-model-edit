#!/usr/bin/env bash
set -euo pipefail

start_port="${1:-}"
if [[ -z "$start_port" ]]; then
  echo "usage: find_free_port.sh <start_port>" >&2
  exit 1
fi

port_in_use() {
  if command -v ss >/dev/null 2>&1; then
    ss -ltn "sport = :$1" | awk 'NR>1 {found=1} END {exit !found}'
  elif command -v lsof >/dev/null 2>&1; then
    lsof -nP -iTCP:"$1" -sTCP:LISTEN >/dev/null 2>&1
  else
    return 1
  fi
}

find_free_port() {
  local candidate="$1"
  for _ in $(seq 0 20); do
    if ! port_in_use "$candidate"; then
      echo "$candidate"
      return 0
    fi
    candidate=$((candidate + 1))
  done
  echo "no free port found (start=$1)" >&2
  return 1
}

find_free_port "$start_port"
