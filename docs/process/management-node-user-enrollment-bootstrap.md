# Management Node User Enrollment Bootstrap

This runbook provides a one-command setup flow for management-node user enrollment.

## Purpose

- Install the latest (or pinned) `bootstrap-management-node.sh` helper.
- Verify script integrity with SHA-256 checksum by default.
- Run enrollment in safe mode by default (dry-run unless `--apply` is passed).

## One-Command Quick Start

Run from any machine with `bash` and `curl`:

```bash
curl -fsSL https://raw.githubusercontent.com/grace-bidos/chem-model-edit/main/scripts/install-management-node.sh \
  | bash -s -- --management-host localhost --user chemops
```

The command above executes only a dry-run plan. No remote changes are applied unless you add `--apply`.

## Version Pinning

Pin to a specific branch, tag, or commit by setting `MANAGEMENT_NODE_VERSION`:

```bash
curl -fsSL https://raw.githubusercontent.com/grace-bidos/chem-model-edit/main/scripts/install-management-node.sh \
  | MANAGEMENT_NODE_VERSION=v0.1.0 bash -s -- --management-host localhost --user chemops
```

You can also set a commit SHA:

```bash
curl -fsSL https://raw.githubusercontent.com/grace-bidos/chem-model-edit/main/scripts/install-management-node.sh \
  | MANAGEMENT_NODE_VERSION=<commit-sha> bash -s -- --management-host localhost --user chemops
```

## Checksum Verification Policy

- Default: checksum verification is enabled.
- Source checksum file: `scripts/bootstrap-management-node.sh.sha256`.
- Disable only for emergency debugging:

```bash
curl -fsSL https://raw.githubusercontent.com/grace-bidos/chem-model-edit/main/scripts/install-management-node.sh \
  | bash -s -- --no-verify-checksum --management-host localhost --user chemops
```

## Apply Mode

Dry-run is default. To apply changes, pass `--apply` through to bootstrap:

```bash
curl -fsSL https://raw.githubusercontent.com/grace-bidos/chem-model-edit/main/scripts/install-management-node.sh \
  | bash -s -- --management-host localhost --user chemops --apply
```

## Local Repository Usage

When working from a checked-out repository:

```bash
scripts/install-management-node.sh --management-host localhost --user chemops
```

## Installer Environment Variables

- `MANAGEMENT_NODE_VERSION` (default: `main`)
- `MANAGEMENT_NODE_REPOSITORY` (default: `grace-bidos/chem-model-edit`)
- `MANAGEMENT_NODE_INSTALL_DIR` (default: `/tmp/management-node-enrollment`)
- `MANAGEMENT_NODE_VERIFY_CHECKSUM` (default: `true`)
