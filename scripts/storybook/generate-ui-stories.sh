#!/usr/bin/env bash
set -euo pipefail

root_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ui_dir="$root_dir/apps/web/src/components/ui"

if [[ ! -d "$ui_dir" ]]; then
  echo "UI directory not found: $ui_dir" >&2
  exit 1
fi

shopt -s nullglob
for file_path in "$ui_dir"/*.tsx; do
  file_name="$(basename "$file_path")"

  if [[ "$file_name" == *.stories.tsx ]]; then
    continue
  fi
  if [[ "$file_name" == *.test.tsx ]]; then
    continue
  fi
  if [[ "$file_name" == *.a11y.test.tsx ]]; then
    continue
  fi
  if [[ "$file_name" == *.fastcheck.test.tsx ]]; then
    continue
  fi

  base_name="${file_name%.tsx}"
  story_path="$ui_dir/$base_name.stories.tsx"

  if [[ -f "$story_path" ]]; then
    echo "skip (already exists): $story_path"
    continue
  fi

  export_name="$(echo "$base_name" | sed -E 's/(^|[-_])([a-z])/\U\2/g')"

  (
    cd "$root_dir"
    pnpm exec hygen story new \
      --name "$base_name" \
      --dir "apps/web/src/components/ui" \
      --exportStyle "named" \
      --exportName "$export_name"
  )
done
