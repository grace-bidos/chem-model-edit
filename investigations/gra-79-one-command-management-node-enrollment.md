# GRA-79 - One-command management-node enrollment setup

## Goal

Provide an operator-friendly one-command UX to run management-node user enrollment with safe defaults and checksum protection.

## Scope

- Add installer script: `scripts/install-management-node.sh`.
- Add checksum artifact: `scripts/bootstrap-management-node.sh.sha256`.
- Keep `scripts/bootstrap-management-node.sh` behavior unchanged and safe by default.
- Document usage and guardrails.

## Decisions

### 1) Keep bootstrap as execution source of truth

The installer only downloads, verifies, and executes `bootstrap-management-node.sh`. All enrollment logic remains in bootstrap.

### 2) Checksum verification default ON

Installer downloads both files:

- `bootstrap-management-node.sh`
- `bootstrap-management-node.sh.sha256`

Then it runs `sha256sum -c` unless explicitly disabled with `--no-verify-checksum` or `MANAGEMENT_NODE_VERIFY_CHECKSUM=false`.

### 3) Version pin support

Installer supports pinning via `MANAGEMENT_NODE_VERSION` and resolves files from:

`https://raw.githubusercontent.com/<repo>/<version>/scripts/...`

This works for branch, tag, or commit SHA pins.

### 4) Safe execution preserved

Installer forwards unknown arguments directly to bootstrap. Since bootstrap defaults to dry-run, the one-command flow is safe by default and requires explicit `--apply` to mutate systems.

## Verification Plan

- Syntax validation:
  - `bash -n scripts/install-management-node.sh scripts/bootstrap-management-node.sh`
- Enrollment harness (if feasible):
  - `scripts/verify-management-node-enrollment-docker.sh`

## Risks and Notes

- The installer assumes `curl` or `wget` is available.
- Checksum verification requires `sha256sum`.
- Network access is required for raw GitHub downloads when not running from local repo paths.
