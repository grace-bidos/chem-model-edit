# Management Node Docker Verification Harness

This runbook provides a local Docker harness to verify management-node bootstrap prerequisites, with Ansible available in-container.

## Scope

- Brings up a disposable container that represents a management node.
- Verifies Ansible and Python tooling in that container.
- Executes a minimal local Ansible ping against `localhost`.

## Prerequisites

- Docker Engine with Compose plugin (`docker compose`)
- Repository cloned locally (any branch)

## File Layout

- Compose file: `ops/docker/docker-compose.management-node.yml`
- Docker build context: `ops/docker/management-node`
- Ansible config: `ops/docker/management-node/ansible.cfg`
- Inventory: `ops/docker/management-node/inventory.ini`

## Commands

### Up

```bash
docker compose -f ops/docker/docker-compose.management-node.yml up -d --build
```

### Verify

```bash
docker compose -f ops/docker/docker-compose.management-node.yml ps

docker compose -f ops/docker/docker-compose.management-node.yml exec management-node ansible --version

docker compose -f ops/docker/docker-compose.management-node.yml exec management-node python --version

docker compose -f ops/docker/docker-compose.management-node.yml exec management-node \
  ansible management -m ping
```

Expected behavior:

- Service `management-node` is `running`
- `ansible --version` prints a valid Ansible release
- `python --version` reports Python 3.12+
- `ansible management -m ping` returns `pong` for `localhost`

### Down

```bash
docker compose -f ops/docker/docker-compose.management-node.yml down
```

### Reset

```bash
docker compose -f ops/docker/docker-compose.management-node.yml down --volumes --remove-orphans
```

If you need a clean image rebuild after reset:

```bash
docker compose -f ops/docker/docker-compose.management-node.yml build --no-cache
```
