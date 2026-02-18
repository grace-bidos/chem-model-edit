#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
cd "$repo_root"

git config core.hooksPath .githooks
echo "Installed local git hooks: core.hooksPath=.githooks"
echo "To bypass once: SKIP_PREPUSH_STRICT=1 git push"
