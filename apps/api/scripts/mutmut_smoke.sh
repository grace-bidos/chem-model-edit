#!/usr/bin/env bash
set -euo pipefail

status=0
uv run mutmut run || status=$?

if [[ -d "mutants" ]]; then
  ts="$(date +%Y%m%d%H%M%S)"
  target=".mutants_local_tmp.${ts}"
  mv mutants "${target}"
  echo "mutmut artifacts moved to ${target}"
fi

exit "${status}"
