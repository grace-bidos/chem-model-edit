#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
COMPOSE_FILE="$ROOT_DIR/ops/ansible/management-node-enrollment/docker-harness/docker-compose.yml"
LOG_DIR="$ROOT_DIR/.just-runtime/management-node-enrollment/verification-logs"
CONTAINER_SERVICE="management-node"

BOOTSTRAP_ARGS="--management-host localhost --user chemops --ssh-public-key 'ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIGRA82EXAMPLE management-admin@gra82'"

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

run_compose() {
  docker compose -f "$COMPOSE_FILE" "$@"
}

run_in_container() {
  run_compose exec -T "$CONTAINER_SERVICE" bash -lc "$1"
}

timestamp() {
  date '+%Y%m%d-%H%M%S'
}

mkdir -p "$LOG_DIR"
STAMP="$(timestamp)"
DRY_LOG="$LOG_DIR/${STAMP}-dry-run.log"
APPLY_LOG="$LOG_DIR/${STAMP}-apply.log"
IDEMPOTENCY_LOG="$LOG_DIR/${STAMP}-idempotency.log"

have_cmd docker || { echo "error: docker command not found" >&2; exit 1; }
run_compose version >/dev/null

echo "[1/6] Reset docker harness"
run_compose down --volumes --remove-orphans >/dev/null 2>&1 || true

echo "[2/6] Build and start docker harness"
run_compose up -d --build

echo "[3/6] Tooling sanity checks"
run_in_container "ansible --version"
run_in_container "python --version"
run_in_container "cd /workspace && bash -n scripts/bootstrap-management-node.sh"

echo "[4/6] Dry-run enrollment"
run_in_container "cd /workspace && scripts/bootstrap-management-node.sh ${BOOTSTRAP_ARGS}" | tee "$DRY_LOG"

echo "[5/6] Apply enrollment"
run_in_container "cd /workspace && scripts/bootstrap-management-node.sh --apply ${BOOTSTRAP_ARGS}" | tee "$APPLY_LOG"

echo "[6/6] Idempotency apply"
run_in_container "cd /workspace && scripts/bootstrap-management-node.sh --apply ${BOOTSTRAP_ARGS}" | tee "$IDEMPOTENCY_LOG"

if ! grep -Eq 'changed=0\s+unreachable=0\s+failed=0' "$IDEMPOTENCY_LOG"; then
  echo "error: idempotency check failed (expected changed=0 in second apply)" >&2
  exit 1
fi

echo "verification_logs:"
echo "  dry-run:     $DRY_LOG"
echo "  apply:       $APPLY_LOG"
echo "  idempotency: $IDEMPOTENCY_LOG"

echo "done"
