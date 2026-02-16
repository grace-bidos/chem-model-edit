# AiiDA VM Bootstrap (GRA-88)

This runbook defines the VM bootstrap flow for PostgreSQL + RabbitMQ + AiiDA profile setup.

Scope of this slice (`GRA-88`):

- Add an idempotent bootstrap script for Ubuntu-like VMs.
- Keep dry-run as the default safety mode.
- Provide explicit apply flow and sanity checks.

## Files

- Bootstrap entrypoint: `scripts/aiida-vm-bootstrap.sh`
- Bootstrap implementation: `ops/aiida-vm/bootstrap-aiida-vm.sh`
- Environment template: `ops/aiida-vm/aiida-vm.env.example`

## Prerequisites

- Ubuntu/Debian VM with `sudo` access
- `uv` installed
- AiiDA dependency group available in API project (`apps/api`)

## 1) Prepare Environment

Copy and edit env file:

```bash
cp ops/aiida-vm/aiida-vm.env.example ops/aiida-vm/aiida-vm.env
```

Minimum values to review:

- `AIIDA_PGPASSWORD`
- `AIIDA_BROKER_PASSWORD`
- `AIIDA_PROFILE`
- `AIIDA_CONFIG_DIR` and `AIIDA_REPOSITORY_DIR`

## 2) Dry-run (required first)

Dry-run prints all commands without applying changes:

```bash
scripts/aiida-vm-bootstrap.sh \
  --copy-env \
  --env-file ops/aiida-vm/aiida-vm.env \
  --init-profile \
  --sanity-check
```

Dry-run covers:

- package/service operations for PostgreSQL and RabbitMQ
- PostgreSQL role/database creation commands
- RabbitMQ user/vhost/permission commands
- AiiDA profile setup command
- sanity check commands (`verdi profile list`, `verdi status`)

## 3) Apply

Execute bootstrap on VM:

```bash
scripts/aiida-vm-bootstrap.sh \
  --apply \
  --env-file ops/aiida-vm/aiida-vm.env \
  --init-profile \
  --sanity-check
```

Behavior in apply mode:

- runs `apt-get update` and installs `postgresql` + `rabbitmq-server` (default)
- enables/starts both services
- waits for TCP readiness (`5432`, `5672`)
- ensures PostgreSQL role/database exist (idempotent)
- ensures RabbitMQ vhost/user/permissions exist (idempotent)
- creates AiiDA profile only when missing
- runs sanity checks

## 4) Manual Sanity Checks

If you need to re-run checks separately:

```bash
AIIDA_PATH=.just-runtime/aiida-vm/config \
uv run --project apps/api --group aiida verdi profile list

AIIDA_PATH=.just-runtime/aiida-vm/config \
uv run --project apps/api --group aiida verdi status
```

## Safety Model

- Dry-run by default; destructive operations are not used.
- Existing DB/broker resources are reused where possible.
- Re-running apply is safe and updates credentials/permissions idempotently.
