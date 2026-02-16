#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  scripts/runner/run_ephemeral_runner.sh \
    --jit-config /tmp/jit-config.json \
    --runner-home /opt/actions-runner \
    [--runner-version 2.327.1]

Required:
  --jit-config <json file from request_jit_config.sh>
  --runner-home <directory to store runner binaries>

Notes:
  - Runs exactly one job when JIT config is valid.
  - On first install, you may need OS dependencies:
      sudo ./bin/installdependencies.sh
EOF
}

jit_config_file=""
runner_home=""
runner_version="2.327.1"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --jit-config) jit_config_file="${2:-}"; shift 2 ;;
    --runner-home) runner_home="${2:-}"; shift 2 ;;
    --runner-version) runner_version="${2:-}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${jit_config_file}" || -z "${runner_home}" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if [[ ! -f "${jit_config_file}" ]]; then
  echo "JIT config file not found: ${jit_config_file}" >&2
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "jq is required." >&2
  exit 1
fi

encoded_jit_config="$(jq -r '.encoded_jit_config // empty' "${jit_config_file}")"
if [[ -z "${encoded_jit_config}" ]]; then
  echo "encoded_jit_config not found in ${jit_config_file}" >&2
  exit 1
fi

mkdir -p "${runner_home}"
cd "${runner_home}"

if [[ ! -x "./run.sh" ]]; then
  archive="actions-runner-linux-x64-${runner_version}.tar.gz"
  url="https://github.com/actions/runner/releases/download/v${runner_version}/${archive}"
  curl -fsSL "${url}" -o "${archive}"
  tar xzf "${archive}"
  rm -f "${archive}"
fi

./run.sh --jitconfig "${encoded_jit_config}"
