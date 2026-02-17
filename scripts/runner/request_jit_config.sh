#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/runner/request_jit_config.sh \
    --scope repo \
    --owner <org-or-user> \
    --repo <repo-name> \
    --runner-group-id <id> \
    --labels "self-hosted,linux,x64,chem-trusted-pr" \
    --name-prefix chem-jit \
    --token-source env \
    --token-env-var GH_TOKEN \
    --out /tmp/jit-config.json

Required:
  --scope repo|org
  --owner <org-or-user>
  --runner-group-id <id>
  --labels <comma-separated labels>
  --out <output json path>
  --token-source gh|env|file        Optional. Default: gh
  --token-env-var <name>            Optional. Default: GH_TOKEN (used for env mode)
  --token-file <path>               Required only for file mode

When --scope repo:
  --repo <repo-name> is required.

Notes:
  - Requires either authenticated `gh` OR explicit token source.
  - Output file contains sensitive encoded JIT config. Treat as secret.
  - Token source can be used with GitHub App installation tokens.
EOF
}

scope=""
owner=""
repo=""
runner_group_id=""
labels_csv=""
name_prefix="chem-jit"
out=""
work_folder="_work"
token_source="gh"
token_file=""
token_env_var="GH_TOKEN"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --scope) scope="${2:-}"; shift 2 ;;
    --owner) owner="${2:-}"; shift 2 ;;
    --repo) repo="${2:-}"; shift 2 ;;
    --runner-group-id) runner_group_id="${2:-}"; shift 2 ;;
    --labels) labels_csv="${2:-}"; shift 2 ;;
    --name-prefix) name_prefix="${2:-}"; shift 2 ;;
    --work-folder) work_folder="${2:-}"; shift 2 ;;
    --out) out="${2:-}"; shift 2 ;;
    --token-source) token_source="${2:-}"; shift 2 ;;
    --token-file) token_file="${2:-}"; shift 2 ;;
    --token-env-var) token_env_var="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$scope" || -z "$owner" || -z "$runner_group_id" || -z "$labels_csv" || -z "$out" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if [[ "$scope" == "repo" && -z "$repo" ]]; then
  echo "--repo is required for repo scope." >&2
  exit 1
fi

if [[ "$scope" != "repo" && "$scope" != "org" ]]; then
  echo "--scope must be repo or org." >&2
  exit 1
fi

if [[ "$token_source" != "gh" && "$token_source" != "env" && "$token_source" != "file" ]]; then
  echo "--token-source must be one of: gh, env, file." >&2
  exit 1
fi

if ! command -v gh >/dev/null 2>&1; then
  echo "gh CLI is required." >&2
  exit 1
fi

gh_token=""
if [[ "$token_source" == "gh" ]]; then
  if ! gh auth status -h github.com >/dev/null 2>&1; then
    echo "gh is not authenticated. Run: gh auth login or use --token-source env|file." >&2
    exit 1
  fi
elif [[ "$token_source" == "env" ]]; then
  gh_token="${!token_env_var:-}"
  if [[ -z "$gh_token" ]]; then
    echo "Token env var is empty: ${token_env_var}" >&2
    exit 1
  fi
elif [[ "$token_source" == "file" ]]; then
  if [[ -z "$token_file" ]]; then
    echo "--token-file is required when --token-source file." >&2
    exit 1
  fi
  if [[ ! -f "$token_file" ]]; then
    echo "Token file not found: ${token_file}" >&2
    exit 1
  fi
  gh_token="$(<"$token_file")"
  gh_token="${gh_token//$'\r'/}"
  gh_token="${gh_token//$'\n'/}"
  if [[ -z "$gh_token" ]]; then
    echo "Token file is empty: ${token_file}" >&2
    exit 1
  fi
fi

timestamp="$(date +%Y%m%d-%H%M%S)"
runner_name="${name_prefix}-${timestamp}"

# Convert CSV labels into JSON array.
labels_json="$(printf '%s' "$labels_csv" | awk -F',' '
  BEGIN { printf "[" }
  {
    for (i=1; i<=NF; i++) {
      gsub(/^[ \t]+|[ \t]+$/, "", $i)
      if ($i != "") {
        if (c++) printf ","
        gsub(/"/, "\\\"", $i)
        printf "\"%s\"", $i
      }
    }
  }
  END { printf "]" }
')"

mkdir -p "$(dirname "$out")"

if [[ "$scope" == "repo" ]]; then
  endpoint="repos/${owner}/${repo}/actions/runners/generate-jitconfig"
else
  endpoint="orgs/${owner}/actions/runners/generate-jitconfig"
fi

gh_api_cmd=(gh api \
  --method POST \
  -H "Accept: application/vnd.github+json" \
  "${endpoint}" \
  -f name="${runner_name}" \
  -F runner_group_id="${runner_group_id}" \
  -f work_folder="${work_folder}" \
  -f labels="${labels_json}")

if [[ "$token_source" == "gh" ]]; then
  "${gh_api_cmd[@]}" > "${out}"
else
  GH_TOKEN="$gh_token" "${gh_api_cmd[@]}" > "${out}"
fi

chmod 600 "${out}" || true
echo "Wrote JIT config to: ${out}"
echo "Runner name: ${runner_name}"
