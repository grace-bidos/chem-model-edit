#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  scripts/runner/check_local_runner_health.sh --owner <org-or-user> --repo <repo-name> [options]

Required:
  --owner <org-or-user>
  --repo <repo-name>

Options:
  --labels <comma-separated labels>   Filter GitHub runners that include all labels.
  --strict-gh                         Treat GitHub online runner without local active unit as degraded.
  -h, --help                          Show this help.

Behavior:
  - Reports local `actions.runner.*` systemd service states.
  - Reports GitHub runner online/offline status for the target repository.
  - Flags local-vs-GitHub status mismatch.
  - Exit code: 0 healthy, non-zero degraded.

Dependencies:
  - systemctl
  - gh (authenticated with repo/admin access for Actions runners)
  - jq
USAGE
}

owner=""
repo=""
labels_csv=""
strict_gh=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --owner) owner="${2:-}"; shift 2 ;;
    --repo) repo="${2:-}"; shift 2 ;;
    --labels) labels_csv="${2:-}"; shift 2 ;;
    --strict-gh) strict_gh=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "$owner" || -z "$repo" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

for tool in systemctl gh jq; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "Missing required tool: $tool" >&2
    exit 2
  fi
done

mapfile -t units < <(systemctl list-units 'actions.runner.*.service' --all --plain --no-legend --no-pager | awk '{print $1}')

declare -A local_state=()
declare -A local_active=()

echo "== Local systemd: actions.runner.* =="
if [[ ${#units[@]} -eq 0 ]]; then
  echo "(no matching units found)"
else
  printf '%-70s %-10s %-10s\n' "UNIT" "ACTIVE" "SUB"
  for unit in "${units[@]}"; do
    mapfile -t states < <(systemctl show "$unit" --property=ActiveState,SubState --value --no-pager)
    active="${states[0]:-unknown}"
    sub="${states[1]:-unknown}"
    local_state["$unit"]="$active/$sub"
    if [[ "$active" == "active" ]]; then
      local_active["$unit"]=1
    fi
    printf '%-70s %-10s %-10s\n' "$unit" "$active" "$sub"
  done
fi

labels_json='[]'
if [[ -n "$labels_csv" ]]; then
  labels_json="$(printf '%s' "$labels_csv" | awk -F',' '
    BEGIN { printf "[" }
    {
      for (i=1; i<=NF; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", $i)
        if ($i != "") {
          if (c++) printf ","
          gsub(/\"/, "\\\"", $i)
          printf "\"%s\"", $i
        }
      }
    }
    END { printf "]" }
  ')"
fi

gh_json="$(gh api "repos/${owner}/${repo}/actions/runners?per_page=100")"

echo
echo "== GitHub runners (repo: ${owner}/${repo}) =="

mapfile -t gh_rows < <(jq -r --argjson required "$labels_json" '
  .runners[]
  | select(
      ($required | length) == 0
      or ((.labels | map(.name)) as $have | ($required - $have | length) == 0)
    )
  | [.name, .status, .busy, (.labels | map(.name) | join(","))]
  | @tsv
' <<<"$gh_json")

if [[ ${#gh_rows[@]} -eq 0 ]]; then
  echo "(no runners found for filter)"
else
  printf '%-40s %-10s %-6s %s\n' "RUNNER" "STATUS" "BUSY" "LABELS"
  for row in "${gh_rows[@]}"; do
    IFS=$'\t' read -r name status busy labels <<<"$row"
    printf '%-40s %-10s %-6s %s\n' "$name" "$status" "$busy" "$labels"
  done
fi

degraded=0
echo
echo "== Mismatch checks =="

# 1) Local active service should not map to offline/missing GitHub runner.
if [[ ${#units[@]} -gt 0 ]]; then
  for unit in "${units[@]}"; do
    if [[ -z "${local_active[$unit]:-}" ]]; then
      continue
    fi

    # Unit format is usually: actions.runner.<scope>.<runner-name>.service
    runner_name="${unit#actions.runner.}"
    runner_name="${runner_name#*.}"
    runner_name="${runner_name%.service}"

    gh_status="$(jq -r --arg n "$runner_name" '.runners[] | select(.name == $n) | .status' <<<"$gh_json" | head -n1)"

    if [[ -z "$gh_status" ]]; then
      echo "DEGRADED: local active unit '${unit}' has no matching GitHub runner '${runner_name}'."
      degraded=1
      continue
    fi

    if [[ "$gh_status" != "online" ]]; then
      echo "DEGRADED: local active unit '${unit}' maps to GitHub status '${gh_status}'."
      degraded=1
    fi
  done
fi

# 2) Optional strict mode: any GitHub online runner should have local active service.
if [[ "$strict_gh" -eq 1 ]]; then
  has_local_active_runner() {
    local target="$1"
    local unit runner_name
    for unit in "${units[@]}"; do
      if [[ -z "${local_active[$unit]:-}" ]]; then
        continue
      fi
      runner_name="${unit#actions.runner.}"
      runner_name="${runner_name#*.}"
      runner_name="${runner_name%.service}"
      if [[ "$runner_name" == "$target" ]]; then
        return 0
      fi
    done
    return 1
  }

  while IFS= read -r gh_online; do
    [[ -z "$gh_online" ]] && continue
    if ! has_local_active_runner "$gh_online"; then
      echo "DEGRADED: GitHub runner '${gh_online}' is online but no local active actions.runner.* service was found."
      degraded=1
    fi
  done < <(jq -r '.runners[] | select(.status == "online") | .name' <<<"$gh_json")
fi

if [[ "$degraded" -eq 0 ]]; then
  echo "HEALTHY: local and GitHub runner states are consistent."
  exit 0
fi

echo "DEGRADED: mismatch detected."
exit 1
