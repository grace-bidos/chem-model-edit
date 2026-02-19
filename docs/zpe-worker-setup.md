# Compute Node Setup (runtime-only)

This path is Redis-free.
Compute node registration uses `/api/runtime/*` and Convex-backed node metadata.

## Prerequisites

- Control-plane API is reachable (`/api/runtime/nodes/*`)
- Clerk sign-in available in Web SaaS
- Ubuntu compute node has `curl` and `jq`

## Operator Flow

1. In Web SaaS, sign in and generate a compute-node setup command.
2. Copy the command shown in UI.
3. Run it on the compute node terminal.

Generated command shape:

```bash
curl -fsSL <API_BASE>/api/runtime/nodes/install.sh | bash -s -- \
  --api-base <API_BASE> \
  --join-token <JOIN_TOKEN> \
  --queue-name <QUEUE_NAME>
```

The installer calls:

- `POST /api/runtime/nodes/register`

and creates a runtime target for the issuing user.

## Notes

- Join token is short-lived and single-use.
- Active target selection is managed via `/api/runtime/targets`.
