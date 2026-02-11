#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <title> <body-file> [label ...]" >&2
  exit 1
fi

title="$1"
body_file="$2"
shift 2

if [[ ! -f "$body_file" ]]; then
  echo "error: body file not found: $body_file" >&2
  exit 1
fi

args=(issue create --title "$title" --body-file "$body_file")
for label in "$@"; do
  args+=(--label "$label")
done

gh "${args[@]}"
