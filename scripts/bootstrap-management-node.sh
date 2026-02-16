#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_ANSIBLE_DIR="$ROOT_DIR/ops/ansible/management-node-enrollment"
DEFAULT_RUNTIME_DIR="$ROOT_DIR/.just-runtime/management-node-enrollment"

DRY_RUN=1
PLAN_ONLY=0
MANAGEMENT_HOST="localhost"
ENROLL_USER="${USER:-chemops}"
SSH_PUBLIC_KEY=""
SUDO_NOPASSWD="true"
ANSIBLE_DIR="$DEFAULT_ANSIBLE_DIR"
INVENTORY_FILE=""
VARS_FILE=""

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Bootstrap management-node user enrollment through Ansible.
Safe by default: dry-run (ansible --check --diff) unless --apply is set.

Options:
  --apply                     Execute changes on target host
  --plan-only                 Generate files and print command without running ansible
  --management-host <host>    Target management node host (default: localhost)
  --user <name>               User account to enroll (default: current shell user)
  --ssh-public-key <value>    Public key literal or @/path/to/key.pub
  --sudo-nopasswd <bool>      true|false (default: true)
  --ansible-dir <path>        Ansible project root (default: ops/ansible/management-node-enrollment)
  --inventory-file <path>     Existing inventory file (skip generated inventory)
  --vars-file <path>          Existing vars file (skip generated vars)
  -h, --help                  Show this help
USAGE
}

log() {
  printf '%s\n' "$*"
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

parse_bool() {
  case "$1" in
    true|false) printf '%s' "$1" ;;
    *) die "invalid boolean: $1 (expected true or false)" ;;
  esac
}

resolve_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

resolve_public_key() {
  local value="$1"

  if [[ -n "$value" ]]; then
    if [[ "$value" == @* ]]; then
      local key_path="${value#@}"
      [[ -f "$key_path" ]] || die "public key file not found: $key_path"
      cat "$key_path"
      return 0
    fi
    printf '%s\n' "$value"
    return 0
  fi

  local candidates=(
    "$HOME/.ssh/id_ed25519.pub"
    "$HOME/.ssh/id_rsa.pub"
  )

  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -f "$candidate" ]]; then
      cat "$candidate"
      return 0
    fi
  done

  die "no SSH public key found; pass --ssh-public-key"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)
      DRY_RUN=0
      ;;
    --plan-only)
      PLAN_ONLY=1
      ;;
    --management-host)
      shift
      [[ $# -gt 0 ]] || die "--management-host requires a value"
      MANAGEMENT_HOST="$1"
      ;;
    --user)
      shift
      [[ $# -gt 0 ]] || die "--user requires a value"
      ENROLL_USER="$1"
      ;;
    --ssh-public-key)
      shift
      [[ $# -gt 0 ]] || die "--ssh-public-key requires a value"
      SSH_PUBLIC_KEY="$1"
      ;;
    --sudo-nopasswd)
      shift
      [[ $# -gt 0 ]] || die "--sudo-nopasswd requires a value"
      SUDO_NOPASSWD="$(parse_bool "$1")"
      ;;
    --ansible-dir)
      shift
      [[ $# -gt 0 ]] || die "--ansible-dir requires a value"
      ANSIBLE_DIR="$1"
      ;;
    --inventory-file)
      shift
      [[ $# -gt 0 ]] || die "--inventory-file requires a value"
      INVENTORY_FILE="$1"
      ;;
    --vars-file)
      shift
      [[ $# -gt 0 ]] || die "--vars-file requires a value"
      VARS_FILE="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      die "unknown option: $1"
      ;;
  esac
  shift
done

ANSIBLE_DIR="$(resolve_abs_path "$ANSIBLE_DIR")"
[[ -d "$ANSIBLE_DIR" ]] || die "ansible dir not found: $ANSIBLE_DIR"

if [[ -n "$INVENTORY_FILE" ]]; then
  INVENTORY_FILE="$(resolve_abs_path "$INVENTORY_FILE")"
  [[ -f "$INVENTORY_FILE" ]] || die "inventory file not found: $INVENTORY_FILE"
fi

if [[ -n "$VARS_FILE" ]]; then
  VARS_FILE="$(resolve_abs_path "$VARS_FILE")"
  [[ -f "$VARS_FILE" ]] || die "vars file not found: $VARS_FILE"
fi

PUBLIC_KEY_CONTENT="$(resolve_public_key "$SSH_PUBLIC_KEY")"
mkdir -p "$DEFAULT_RUNTIME_DIR"

if [[ -z "$INVENTORY_FILE" ]]; then
  INVENTORY_FILE="$DEFAULT_RUNTIME_DIR/inventory.ini"
  INVENTORY_HOST_LINE="$MANAGEMENT_HOST ansible_user=root"
  if [[ "$MANAGEMENT_HOST" == "localhost" || "$MANAGEMENT_HOST" == "127.0.0.1" || "$MANAGEMENT_HOST" == "::1" ]]; then
    INVENTORY_HOST_LINE="$MANAGEMENT_HOST ansible_connection=local ansible_user=root"
  fi
  cat > "$INVENTORY_FILE" <<INVENTORY
[management_nodes]
$INVENTORY_HOST_LINE
INVENTORY
fi

if [[ -z "$VARS_FILE" ]]; then
  VARS_FILE="$DEFAULT_RUNTIME_DIR/enrollment.vars.yml"
  cat > "$VARS_FILE" <<VARS
---
management_enrollment_users:
  - name: $ENROLL_USER
    groups:
      - sudo
    shell: /bin/bash
    ssh_public_key: "$PUBLIC_KEY_CONTENT"
    sudo_nopasswd: $SUDO_NOPASSWD
VARS
fi

PLAYBOOK_PATH="$ANSIBLE_DIR/playbooks/enroll-users.yml"
[[ -f "$PLAYBOOK_PATH" ]] || die "playbook not found: $PLAYBOOK_PATH"

CMD=(env "ANSIBLE_CONFIG=$ANSIBLE_DIR/ansible.cfg" ansible-playbook -i "$INVENTORY_FILE" "$PLAYBOOK_PATH" -e "@$VARS_FILE")
if [[ "$DRY_RUN" -eq 1 ]]; then
  CMD+=(--check --diff)
fi

log "Management-node enrollment bootstrap"
log "  Mode:             $([[ "$DRY_RUN" -eq 1 ]] && echo dry-run || echo apply)"
log "  Plan only:        $([[ "$PLAN_ONLY" -eq 1 ]] && echo yes || echo no)"
log "  Management host:  $MANAGEMENT_HOST"
log "  User:             $ENROLL_USER"
log "  Inventory file:   $INVENTORY_FILE"
log "  Vars file:        $VARS_FILE"
log "  Ansible dir:      $ANSIBLE_DIR"

printf 'Command:'
printf ' %q' "${CMD[@]}"
printf '\n'

if [[ "$PLAN_ONLY" -eq 1 ]]; then
  exit 0
fi

command -v ansible-playbook >/dev/null 2>&1 || die "ansible-playbook not found. Install Ansible first or pass --plan-only"

"${CMD[@]}"
