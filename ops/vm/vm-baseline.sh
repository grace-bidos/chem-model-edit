#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PS_SCRIPT_WIN="$(wslpath -w "$SCRIPT_DIR/hyperv-baseline.ps1")"

if ! command -v powershell.exe >/dev/null 2>&1; then
  echo "powershell.exe not found. Run this script from WSL on Windows." >&2
  exit 1
fi

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$PS_SCRIPT_WIN" "$@"
