#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_ANSIBLE_DIR="$ROOT_DIR/ops/ansible/slurm-byo/examples"
DEFAULT_RUNTIME_DIR="$ROOT_DIR/.just-runtime/slurm-byo-ansible-wrapper"

DRY_RUN=1
PLAN_ONLY=0
CLUSTER_NAME="chem-byo"
CONTROLLER_HOST=""
SSH_USER="root"
SSH_PRIVATE_KEY=""
SLURM_PARTITION="compute"
SLURM_ACCOUNT=""
SLURM_QOS=""
ANSIBLE_DIR="$DEFAULT_ANSIBLE_DIR"
PLAYBOOK_RELATIVE="playbooks/site.yml"
INVENTORY_OUT=""
VARS_OUT=""
WORKER_HOSTS=()

usage() {
  cat <<USAGE
Usage: $(basename "$0") [options]

Generate inventory/vars for a minimal Slurm BYO control+worker setup,
then run ansible-playbook in check-mode by default.

Safe defaults:
  - check-mode run (ansible --check --diff)
  - no apply unless --apply

Options:
  --apply                      Execute ansible without --check/--diff
  --plan-only                  Generate files and print command only
  --cluster-name <name>        Slurm ClusterName (default: chem-byo)
  --controller-host <host>     Controller host or IP (required)
  --worker-host <host>         Worker host or IP (repeatable, required >=1)
  --worker-hosts <csv>         Comma-separated worker hosts
  --ssh-user <name>            SSH user for all hosts (default: root)
  --ssh-private-key <path>     SSH private key for ansible_ssh_private_key_file
  --partition <name>           Partition name (default: compute)
  --account <name>             Optional AccountingStorage account hint
  --qos <name>                 Optional QoS hint
  --ansible-dir <path>         Ansible directory (default: ops/ansible/slurm-byo/examples)
  --playbook <path>            Playbook path relative to ansible-dir (default: playbooks/site.yml)
  --inventory-out <path>       Generated inventory path
  --vars-out <path>            Generated vars path
  -h, --help                   Show this help
USAGE
}

log() {
  printf '%s\n' "$*"
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

resolve_abs_path() {
  local value="$1"
  if [[ "$value" = /* ]]; then
    printf '%s\n' "$value"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$value"
  fi
}

append_worker_hosts_csv() {
  local csv="$1"
  local item
  IFS=',' read -r -a items <<<"$csv"
  for item in "${items[@]}"; do
    item="${item## }"
    item="${item%% }"
    if [[ -n "$item" ]]; then
      WORKER_HOSTS+=("$item")
    fi
  done
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --apply)
      DRY_RUN=0
      ;;
    --plan-only)
      PLAN_ONLY=1
      ;;
    --cluster-name)
      shift
      [[ $# -gt 0 ]] || die "--cluster-name requires a value"
      CLUSTER_NAME="$1"
      ;;
    --controller-host)
      shift
      [[ $# -gt 0 ]] || die "--controller-host requires a value"
      CONTROLLER_HOST="$1"
      ;;
    --worker-host)
      shift
      [[ $# -gt 0 ]] || die "--worker-host requires a value"
      WORKER_HOSTS+=("$1")
      ;;
    --worker-hosts)
      shift
      [[ $# -gt 0 ]] || die "--worker-hosts requires a value"
      append_worker_hosts_csv "$1"
      ;;
    --ssh-user)
      shift
      [[ $# -gt 0 ]] || die "--ssh-user requires a value"
      SSH_USER="$1"
      ;;
    --ssh-private-key)
      shift
      [[ $# -gt 0 ]] || die "--ssh-private-key requires a value"
      SSH_PRIVATE_KEY="$1"
      ;;
    --partition)
      shift
      [[ $# -gt 0 ]] || die "--partition requires a value"
      SLURM_PARTITION="$1"
      ;;
    --account)
      shift
      [[ $# -gt 0 ]] || die "--account requires a value"
      SLURM_ACCOUNT="$1"
      ;;
    --qos)
      shift
      [[ $# -gt 0 ]] || die "--qos requires a value"
      SLURM_QOS="$1"
      ;;
    --ansible-dir)
      shift
      [[ $# -gt 0 ]] || die "--ansible-dir requires a value"
      ANSIBLE_DIR="$1"
      ;;
    --playbook)
      shift
      [[ $# -gt 0 ]] || die "--playbook requires a value"
      PLAYBOOK_RELATIVE="$1"
      ;;
    --inventory-out)
      shift
      [[ $# -gt 0 ]] || die "--inventory-out requires a value"
      INVENTORY_OUT="$1"
      ;;
    --vars-out)
      shift
      [[ $# -gt 0 ]] || die "--vars-out requires a value"
      VARS_OUT="$1"
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

[[ -n "$CONTROLLER_HOST" ]] || die "--controller-host is required"
[[ ${#WORKER_HOSTS[@]} -gt 0 ]] || die "at least one --worker-host (or --worker-hosts) is required"
[[ -n "$CLUSTER_NAME" ]] || die "--cluster-name cannot be empty"
[[ -n "$SLURM_PARTITION" ]] || die "--partition cannot be empty"
[[ -n "$SSH_USER" ]] || die "--ssh-user cannot be empty"

ANSIBLE_DIR="$(resolve_abs_path "$ANSIBLE_DIR")"
[[ -d "$ANSIBLE_DIR" ]] || die "ansible dir not found: $ANSIBLE_DIR"

if [[ -n "$SSH_PRIVATE_KEY" ]]; then
  SSH_PRIVATE_KEY="$(resolve_abs_path "$SSH_PRIVATE_KEY")"
  [[ -f "$SSH_PRIVATE_KEY" ]] || die "ssh private key not found: $SSH_PRIVATE_KEY"
fi

mkdir -p "$DEFAULT_RUNTIME_DIR"

if [[ -z "$INVENTORY_OUT" ]]; then
  INVENTORY_OUT="$DEFAULT_RUNTIME_DIR/inventory.ini"
else
  INVENTORY_OUT="$(resolve_abs_path "$INVENTORY_OUT")"
fi

if [[ -z "$VARS_OUT" ]]; then
  VARS_OUT="$DEFAULT_RUNTIME_DIR/group_vars.all.yml"
else
  VARS_OUT="$(resolve_abs_path "$VARS_OUT")"
fi

mkdir -p "$(dirname "$INVENTORY_OUT")" "$(dirname "$VARS_OUT")"

ssh_private_key_kv=""
if [[ -n "$SSH_PRIVATE_KEY" ]]; then
  ssh_private_key_kv=" ansible_ssh_private_key_file=$SSH_PRIVATE_KEY"
fi

{
  printf '[slurmservers]\n'
  printf 'controller ansible_host=%s ansible_user=%s%s\n' "$CONTROLLER_HOST" "$SSH_USER" "$ssh_private_key_kv"
  printf '\n[slurmexechosts]\n'
  worker_index=1
  for worker in "${WORKER_HOSTS[@]}"; do
    printf 'worker%d ansible_host=%s ansible_user=%s%s\n' "$worker_index" "$worker" "$SSH_USER" "$ssh_private_key_kv"
    worker_index=$((worker_index + 1))
  done
  printf '\n[slurm_cluster:children]\n'
  printf 'slurmservers\n'
  printf 'slurmexechosts\n'
} >"$INVENTORY_OUT"

node_lines=""
nodes_csv=""
worker_index=1
for worker in "${WORKER_HOSTS[@]}"; do
  alias="worker${worker_index}"
  node_lines+="  - name: $alias"$'\n'
  node_lines+="    State: UNKNOWN"$'\n'
  if [[ -n "$nodes_csv" ]]; then
    nodes_csv+=","
  fi
  nodes_csv+="$alias"
  worker_index=$((worker_index + 1))
done

{
  printf -- '---\n'
  printf 'slurm_roles:\n'
  printf '  - controller\n'
  printf 'slurm_config:\n'
  printf '  ClusterName: %s\n' "$CLUSTER_NAME"
  printf '  SlurmctldHost: controller\n'
  if [[ -n "$SLURM_ACCOUNT" ]]; then
    printf '  AccountingStorageEnforce: associations\n'
  fi
  printf 'slurm_nodes:\n%s' "$node_lines"
  printf 'slurm_partitions:\n'
  printf '  - name: %s\n' "$SLURM_PARTITION"
  printf '    Default: YES\n'
  printf '    Nodes: "%s"\n' "$nodes_csv"
  if [[ -n "$SLURM_ACCOUNT" ]]; then
    printf 'slurm_byo_account: "%s"\n' "$SLURM_ACCOUNT"
  fi
  if [[ -n "$SLURM_QOS" ]]; then
    printf 'slurm_byo_qos: "%s"\n' "$SLURM_QOS"
  fi
  printf 'slurm_byo_wrapper_metadata:\n'
  printf '  generated_by: bootstrap-slurm-byo.sh\n'
  printf '  generated_at_utc: "%s"\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  mode: "%s"\n' "$([[ "$DRY_RUN" -eq 1 ]] && printf dry-run || printf apply)"
} >"$VARS_OUT"

PLAYBOOK_PATH="$ANSIBLE_DIR/$PLAYBOOK_RELATIVE"
[[ -f "$PLAYBOOK_PATH" ]] || die "playbook not found: $PLAYBOOK_PATH"

CMD=(env "ANSIBLE_CONFIG=$ANSIBLE_DIR/ansible.cfg" ansible-playbook -i "$INVENTORY_OUT" "$PLAYBOOK_PATH" -e "@$VARS_OUT")
if [[ "$DRY_RUN" -eq 1 ]]; then
  CMD+=(--check --diff)
fi

log "Slurm BYO Ansible wrapper PoC"
log "  Mode:            $([[ "$DRY_RUN" -eq 1 ]] && echo dry-run || echo apply)"
log "  Plan only:       $([[ "$PLAN_ONLY" -eq 1 ]] && echo yes || echo no)"
log "  Cluster name:    $CLUSTER_NAME"
log "  Controller host: $CONTROLLER_HOST"
log "  Worker count:    ${#WORKER_HOSTS[@]}"
log "  Inventory file:  $INVENTORY_OUT"
log "  Vars file:       $VARS_OUT"
log "  Ansible dir:     $ANSIBLE_DIR"

printf 'Command:'
printf ' %q' "${CMD[@]}"
printf '\n'

if [[ "$PLAN_ONLY" -eq 1 ]]; then
  exit 0
fi

command -v ansible-playbook >/dev/null 2>&1 || die "ansible-playbook not found. Install Ansible first or run with --plan-only"
"${CMD[@]}"
