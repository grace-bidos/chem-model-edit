#!/usr/bin/env bash
set -euo pipefail

MANAGEMENT_NODE_VERSION="${MANAGEMENT_NODE_VERSION:-main}"
MANAGEMENT_NODE_REPOSITORY="${MANAGEMENT_NODE_REPOSITORY:-grace-bidos/chem-model-edit}"
MANAGEMENT_NODE_INSTALL_DIR="${MANAGEMENT_NODE_INSTALL_DIR:-${TMPDIR:-/tmp}/management-node-enrollment}"
MANAGEMENT_NODE_VERIFY_CHECKSUM="${MANAGEMENT_NODE_VERIFY_CHECKSUM:-true}"
DOWNLOAD_ONLY=0
BOOTSTRAP_ARGS=()

usage() {
  cat <<USAGE
Usage: $(basename "$0") [installer-options] [-- bootstrap-options]

Download and run management-node enrollment bootstrap script.
Checksum verification is enabled by default.

Installer options:
  --version <value>          Set MANAGEMENT_NODE_VERSION (default: main)
  --install-dir <path>       Download directory (default: /tmp/management-node-enrollment)
  --verify-checksum          Enable checksum verification (default)
  --no-verify-checksum       Disable checksum verification
  --download-only            Download files without running bootstrap script
  -h, --help                 Show this help

Environment overrides:
  MANAGEMENT_NODE_VERSION
  MANAGEMENT_NODE_REPOSITORY
  MANAGEMENT_NODE_INSTALL_DIR
  MANAGEMENT_NODE_VERIFY_CHECKSUM

Notes:
  - Bootstrap script remains safe by default and requires --apply for changes.
  - Unknown flags are forwarded to bootstrap-management-node.sh.
USAGE
}

die() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

parse_bool() {
  case "$1" in
    true|false) printf '%s\n' "$1" ;;
    *) die "invalid boolean value: $1 (expected true or false)" ;;
  esac
}

have_cmd() {
  command -v "$1" >/dev/null 2>&1
}

download_file() {
  local url="$1"
  local output="$2"

  if have_cmd curl; then
    curl -fsSL "$url" -o "$output"
    return 0
  fi

  if have_cmd wget; then
    wget -qO "$output" "$url"
    return 0
  fi

  die "missing downloader: install curl or wget"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --version)
      shift
      [[ $# -gt 0 ]] || die "--version requires a value"
      MANAGEMENT_NODE_VERSION="$1"
      ;;
    --install-dir)
      shift
      [[ $# -gt 0 ]] || die "--install-dir requires a value"
      MANAGEMENT_NODE_INSTALL_DIR="$1"
      ;;
    --verify-checksum)
      MANAGEMENT_NODE_VERIFY_CHECKSUM="true"
      ;;
    --no-verify-checksum)
      MANAGEMENT_NODE_VERIFY_CHECKSUM="false"
      ;;
    --download-only)
      DOWNLOAD_ONLY=1
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        BOOTSTRAP_ARGS+=("$1")
        shift
      done
      break
      ;;
    *)
      BOOTSTRAP_ARGS+=("$1")
      ;;
  esac
  shift
done

MANAGEMENT_NODE_VERIFY_CHECKSUM="$(parse_bool "$MANAGEMENT_NODE_VERIFY_CHECKSUM")"
BASE_URL="https://raw.githubusercontent.com/${MANAGEMENT_NODE_REPOSITORY}/${MANAGEMENT_NODE_VERSION}/scripts"

mkdir -p "$MANAGEMENT_NODE_INSTALL_DIR"
BOOTSTRAP_PATH="$MANAGEMENT_NODE_INSTALL_DIR/bootstrap-management-node.sh"
CHECKSUM_PATH="$MANAGEMENT_NODE_INSTALL_DIR/bootstrap-management-node.sh.sha256"

printf 'Downloading bootstrap script from %s\n' "$BASE_URL"
download_file "$BASE_URL/bootstrap-management-node.sh" "$BOOTSTRAP_PATH"
download_file "$BASE_URL/bootstrap-management-node.sh.sha256" "$CHECKSUM_PATH"

if [[ "$MANAGEMENT_NODE_VERIFY_CHECKSUM" == "true" ]]; then
  have_cmd sha256sum || die "sha256sum is required for checksum verification"
  (
    cd "$MANAGEMENT_NODE_INSTALL_DIR"
    sha256sum -c "$(basename "$CHECKSUM_PATH")" --status
  ) || die "checksum verification failed"
  printf 'Checksum verification: passed\n'
else
  printf 'Checksum verification: skipped\n'
fi

chmod +x "$BOOTSTRAP_PATH"

if [[ "$DOWNLOAD_ONLY" -eq 1 ]]; then
  printf 'Downloaded files:\n'
  printf '  %s\n' "$BOOTSTRAP_PATH"
  printf '  %s\n' "$CHECKSUM_PATH"
  exit 0
fi

printf 'Running bootstrap script (safe default is dry-run unless --apply is passed).\n'
exec "$BOOTSTRAP_PATH" "${BOOTSTRAP_ARGS[@]}"
