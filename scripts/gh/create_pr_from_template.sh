#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 ]]; then
  echo "usage: $0 <base> <head> <title> <body-file> [--draft]" >&2
  exit 1
fi

base="$1"
head="$2"
title="$3"
body_file="$4"
shift 4

if [[ ! -f "$body_file" ]]; then
  echo "error: body file not found: $body_file" >&2
  exit 1
fi

args=(pr create --base "$base" --head "$head" --title "$title" --body-file "$body_file")
if [[ "${1:-}" == "--draft" ]]; then
  args+=(--draft)
fi

gh "${args[@]}"
