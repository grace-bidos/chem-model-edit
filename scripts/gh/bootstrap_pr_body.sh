#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
usage: scripts/gh/bootstrap_pr_body.sh <LINEAR-ISSUE> <OUTPUT-MD> [options]

Options:
  --type <Ask|Show|Ship>            Default: Ship
  --size <XS|S|M|L>                 Default: S
  --queue <Required|Optional>       Default: Optional
  --stack <Standalone|Base|Depends on #N>
                                    Default: Standalone
  --coderabbit <Required|Optional>  Default: Optional
  --force                           Overwrite output file if it exists

Example:
  scripts/gh/bootstrap_pr_body.sh GRA-60 /tmp/pr_body.md --type Ship --size XS --queue Optional --stack Standalone --coderabbit Optional
USAGE
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 ]]; then
  usage >&2
  exit 1
fi

linear_issue="$1"
output_file="$2"
shift 2

type_value="Ship"
size_value="S"
queue_policy="Optional"
stack_value="Standalone"
coderabbit_policy="Optional"
force_write="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --type)
      type_value="${2:-}"
      shift 2
      ;;
    --size)
      size_value="${2:-}"
      shift 2
      ;;
    --queue)
      queue_policy="${2:-}"
      shift 2
      ;;
    --stack)
      stack_value="${2:-}"
      shift 2
      ;;
    --coderabbit)
      coderabbit_policy="${2:-}"
      shift 2
      ;;
    --force)
      force_write="true"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! "$linear_issue" =~ ^[A-Z]+-[0-9]+$ ]]; then
  echo "error: LINEAR-ISSUE must look like KEY-123 (received: $linear_issue)" >&2
  exit 1
fi

case "$type_value" in
  Ask|Show|Ship) ;;
  *)
    echo "error: --type must be Ask, Show, or Ship" >&2
    exit 1
    ;;
esac

case "$size_value" in
  XS|S|M|L) ;;
  *)
    echo "error: --size must be XS, S, M, or L" >&2
    exit 1
    ;;
esac

case "$queue_policy" in
  Required|Optional) ;;
  *)
    echo "error: --queue must be Required or Optional" >&2
    exit 1
    ;;
esac

case "$coderabbit_policy" in
  Required|Optional) ;;
  *)
    echo "error: --coderabbit must be Required or Optional" >&2
    exit 1
    ;;
esac

if [[ "$stack_value" != "Standalone" && "$stack_value" != "Base" && ! "$stack_value" =~ ^Depends\ on\ \#[0-9]+$ ]]; then
  echo "error: --stack must be Standalone, Base, or 'Depends on #<PR>'" >&2
  exit 1
fi

if [[ -e "$output_file" && "$force_write" != "true" ]]; then
  echo "error: output file already exists: $output_file (use --force to overwrite)" >&2
  exit 1
fi

required_checked=" "
optional_checked=" "
if [[ "$coderabbit_policy" == "Required" ]]; then
  required_checked="x"
else
  optional_checked="x"
fi

mkdir -p "$(dirname "$output_file")"

cat > "$output_file" <<EOF2
# Summary

- TODO

## Work Item Metadata

- Linear Issue: $linear_issue
- Type: $type_value
- Size: $size_value
- Queue Policy: $queue_policy
- Stack: $stack_value

## Linked Issues

- None

## Changes

- TODO

## Validation

- [ ] Local checks passed
- [ ] Added/updated tests where needed

## Temporary Behavior

- [x] None
- [ ] Present (describe clearly below)
- Description:

## Final Behavior

- TODO

## Follow-up Issue/PR (if any)

- [x] None
- [ ] Required (link issue/PR)
- Link:

## CodeRabbit Policy

- [$required_checked] Required (product/API/auth/worker behavior changed)
- [$optional_checked] Optional (process/docs/template/script-only change)
EOF2

echo "wrote $output_file"
