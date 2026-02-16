#!/usr/bin/env bash
set -uo pipefail

SCRIPT_NAME="preflight-worker-single-node.sh"
SCRIPT_VERSION="1"
OUTPUT_JSON=0

CHECK_COUNT=0
PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0
OVERALL_STATUS="pass"

declare -a CHECK_IDS=()
declare -a CHECK_STATUS=()
declare -a CHECK_MESSAGE=()
declare -a CHECK_DETAILS=()

usage() {
  cat <<'USAGE'
Usage: preflight-worker-single-node.sh [--json] [--help]

Checks single-node worker prerequisites for Apptainer + MPI v1.

Options:
  --json    Emit machine-readable JSON output.
  -h, --help
            Show this help message.

Exit codes:
  0  No failing checks (pass/warn only)
  1  One or more failing checks
USAGE
}

have_command() {
  command -v "$1" >/dev/null 2>&1
}

json_escape() {
  local value="$1"
  value=${value//\\/\\\\}
  value=${value//\"/\\\"}
  value=${value//$'\n'/\\n}
  value=${value//$'\r'/\\r}
  value=${value//$'\t'/\\t}
  printf '%s' "$value"
}

add_check() {
  local id="$1"
  local status="$2"
  local message="$3"
  local details="${4:-}"

  CHECK_IDS+=("$id")
  CHECK_STATUS+=("$status")
  CHECK_MESSAGE+=("$message")
  CHECK_DETAILS+=("$details")
  CHECK_COUNT=$((CHECK_COUNT + 1))

  case "$status" in
    pass)
      PASS_COUNT=$((PASS_COUNT + 1))
      ;;
    warn)
      WARN_COUNT=$((WARN_COUNT + 1))
      if [[ "$OVERALL_STATUS" == "pass" ]]; then
        OVERALL_STATUS="warn"
      fi
      ;;
    fail)
      FAIL_COUNT=$((FAIL_COUNT + 1))
      OVERALL_STATUS="fail"
      ;;
    *)
      ;;
  esac
}

read_file_value() {
  local path="$1"
  if [[ -r "$path" ]]; then
    tr -d '[:space:]' < "$path"
    return 0
  fi
  return 1
}

check_kernel_and_namespaces() {
  local kernel_release
  kernel_release="$(uname -r 2>/dev/null || printf 'unknown')"
  add_check \
    "kernel.release" \
    "pass" \
    "Kernel release detected" \
    "kernel_release=${kernel_release}"

  local max_user_ns
  if max_user_ns="$(read_file_value /proc/sys/user/max_user_namespaces)"; then
    if [[ "$max_user_ns" =~ ^[0-9]+$ ]] && (( max_user_ns > 0 )); then
      add_check \
        "kernel.userns.max_user_namespaces" \
        "pass" \
        "User namespaces are enabled" \
        "max_user_namespaces=${max_user_ns}"
    else
      add_check \
        "kernel.userns.max_user_namespaces" \
        "fail" \
        "User namespaces appear disabled" \
        "max_user_namespaces=${max_user_ns}"
    fi
  else
    add_check \
      "kernel.userns.max_user_namespaces" \
      "warn" \
      "Could not read user namespace limit" \
      "path=/proc/sys/user/max_user_namespaces"
  fi

  local unpriv_userns
  if unpriv_userns="$(read_file_value /proc/sys/kernel/unprivileged_userns_clone)"; then
    if [[ "$unpriv_userns" == "1" ]]; then
      add_check \
        "kernel.userns.unprivileged_clone" \
        "pass" \
        "Unprivileged user namespace clone is enabled" \
        "unprivileged_userns_clone=${unpriv_userns}"
    else
      add_check \
        "kernel.userns.unprivileged_clone" \
        "warn" \
        "Unprivileged user namespace clone is disabled" \
        "unprivileged_userns_clone=${unpriv_userns}; rootless Apptainer features may be limited"
    fi
  else
    add_check \
      "kernel.userns.unprivileged_clone" \
      "warn" \
      "Kernel does not expose unprivileged_userns_clone" \
      "path=/proc/sys/kernel/unprivileged_userns_clone"
  fi
}

check_fuse_and_overlay() {
  if [[ -e /dev/fuse ]]; then
    add_check \
      "fuse.device" \
      "pass" \
      "FUSE device exists" \
      "path=/dev/fuse"
  else
    add_check \
      "fuse.device" \
      "fail" \
      "FUSE device is missing" \
      "path=/dev/fuse"
  fi

  if have_command fusermount3; then
    add_check \
      "fuse.fusermount" \
      "pass" \
      "fusermount3 command is available" \
      "command=fusermount3"
  elif have_command fusermount; then
    add_check \
      "fuse.fusermount" \
      "pass" \
      "fusermount command is available" \
      "command=fusermount"
  else
    add_check \
      "fuse.fusermount" \
      "warn" \
      "fusermount command is not available" \
      "expected=fusermount3|fusermount"
  fi

  if [[ -r /proc/filesystems ]] && grep -Eq '(^|[[:space:]])overlay$' /proc/filesystems; then
    add_check \
      "overlay.filesystem" \
      "pass" \
      "overlay filesystem support is present" \
      "source=/proc/filesystems"
  else
    add_check \
      "overlay.filesystem" \
      "warn" \
      "overlay filesystem support not detected" \
      "source=/proc/filesystems; writable overlay may be unavailable"
  fi
}

check_apptainer() {
  if ! have_command apptainer; then
    add_check \
      "apptainer.command" \
      "fail" \
      "Apptainer command is missing" \
      "install apptainer to run containerized worker jobs"
    add_check \
      "apptainer.exec_help" \
      "warn" \
      "Skipped Apptainer behavior probe because command is missing" \
      "dependency=apptainer"
    return
  fi

  local version_output
  version_output="$(apptainer --version 2>/dev/null || true)"
  if [[ -n "$version_output" ]]; then
    add_check \
      "apptainer.command" \
      "pass" \
      "Apptainer command is available" \
      "version=${version_output}"
  else
    add_check \
      "apptainer.command" \
      "warn" \
      "Apptainer command exists but version probe returned no output" \
      "command=apptainer --version"
  fi

  if apptainer exec --help >/dev/null 2>&1; then
    add_check \
      "apptainer.exec_help" \
      "pass" \
      "Apptainer exec subcommand responds" \
      "probe=apptainer exec --help"
  else
    add_check \
      "apptainer.exec_help" \
      "fail" \
      "Apptainer exec subcommand probe failed" \
      "probe=apptainer exec --help"
  fi
}

check_mpi_multinode_assumptions() {
  local warnings=()

  if [[ -n "${SLURM_NNODES:-}" && "${SLURM_NNODES}" =~ ^[0-9]+$ && "${SLURM_NNODES}" -gt 1 ]]; then
    warnings+=("SLURM_NNODES=${SLURM_NNODES}")
  fi

  if [[ -n "${SLURM_JOB_NUM_NODES:-}" && "${SLURM_JOB_NUM_NODES}" =~ ^[0-9]+$ && "${SLURM_JOB_NUM_NODES}" -gt 1 ]]; then
    warnings+=("SLURM_JOB_NUM_NODES=${SLURM_JOB_NUM_NODES}")
  fi

  if [[ -n "${OMPI_MCA_orte_num_nodes:-}" && "${OMPI_MCA_orte_num_nodes}" =~ ^[0-9]+$ && "${OMPI_MCA_orte_num_nodes}" -gt 1 ]]; then
    warnings+=("OMPI_MCA_orte_num_nodes=${OMPI_MCA_orte_num_nodes}")
  fi

  if [[ -n "${HYDRA_NODELIST:-}" ]]; then
    warnings+=("HYDRA_NODELIST is set")
  fi

  if (( ${#warnings[@]} > 0 )); then
    add_check \
      "mpi.multinode_assumption" \
      "warn" \
      "Detected multinode MPI indicators; v1 worker supports single-node execution only" \
      "$(IFS='; '; printf '%s' "${warnings[*]}")"
  else
    add_check \
      "mpi.multinode_assumption" \
      "pass" \
      "No multinode MPI assumptions detected in current environment" \
      "scope=v1 single-node"
  fi

  if have_command mpirun; then
    add_check \
      "mpi.mpirun.command" \
      "pass" \
      "mpirun command is available" \
      "command=mpirun"
  else
    add_check \
      "mpi.mpirun.command" \
      "warn" \
      "mpirun command is not available" \
      "single-node worker can run without MPI launcher depending on workload"
  fi
}

emit_human() {
  printf 'Worker single-node preflight (Apptainer + MPI v1)\n'
  printf '================================================\n'

  local i
  for ((i = 0; i < CHECK_COUNT; i++)); do
    printf '[%s] %s: %s\n' \
      "${CHECK_STATUS[$i]^^}" \
      "${CHECK_IDS[$i]}" \
      "${CHECK_MESSAGE[$i]}"
    if [[ -n "${CHECK_DETAILS[$i]}" ]]; then
      printf '  details: %s\n' "${CHECK_DETAILS[$i]}"
    fi
  done

  printf '\nSummary: status=%s pass=%d warn=%d fail=%d\n' \
    "$OVERALL_STATUS" \
    "$PASS_COUNT" \
    "$WARN_COUNT" \
    "$FAIL_COUNT"

  if [[ "$OVERALL_STATUS" == "fail" ]]; then
    printf 'Result: FAIL (one or more mandatory prerequisites are missing)\n'
  else
    printf 'Result: OK (safe to continue with caveats above, if any)\n'
  fi
}

emit_json() {
  local now_utc
  now_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  printf '{\n'
  printf '  "checker": "%s",\n' "$(json_escape "$SCRIPT_NAME")"
  printf '  "version": "%s",\n' "$(json_escape "$SCRIPT_VERSION")"
  printf '  "timestamp_utc": "%s",\n' "$(json_escape "$now_utc")"
  printf '  "overall_status": "%s",\n' "$(json_escape "$OVERALL_STATUS")"
  printf '  "counts": {"pass": %d, "warn": %d, "fail": %d},\n' "$PASS_COUNT" "$WARN_COUNT" "$FAIL_COUNT"
  printf '  "checks": [\n'

  local i
  for ((i = 0; i < CHECK_COUNT; i++)); do
    printf '    {"id": "%s", "status": "%s", "message": "%s", "details": "%s"}' \
      "$(json_escape "${CHECK_IDS[$i]}")" \
      "$(json_escape "${CHECK_STATUS[$i]}")" \
      "$(json_escape "${CHECK_MESSAGE[$i]}")" \
      "$(json_escape "${CHECK_DETAILS[$i]}")"
    if (( i < CHECK_COUNT - 1 )); then
      printf ','
    fi
    printf '\n'
  done

  printf '  ]\n'
  printf '}\n'
}

parse_args() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --json)
        OUTPUT_JSON=1
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        printf 'error: unknown option: %s\n\n' "$1" >&2
        usage >&2
        exit 1
        ;;
    esac
    shift
  done
}

main() {
  parse_args "$@"

  check_kernel_and_namespaces
  check_fuse_and_overlay
  check_apptainer
  check_mpi_multinode_assumptions

  if (( OUTPUT_JSON == 1 )); then
    emit_json
  else
    emit_human
  fi

  if [[ "$OVERALL_STATUS" == "fail" ]]; then
    exit 1
  fi
  exit 0
}

main "$@"
