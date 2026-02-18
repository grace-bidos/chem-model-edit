# Convex MCP Authentication Troubleshooting Log (Codex + CLI)

## Document intent

This is an operator-focused troubleshooting record for a real Convex MCP authentication failure in Codex.
It captures the exact symptom pattern, investigation flow, failed hypotheses, and the final stable workaround.

The writing style is intentionally detailed so this can be reused as a technical blog draft.

## Scope and context

- Project: `chem-model-edit`
- Host: local Linux/WSL environment
- Agent runtime: Codex with MCP servers configured in `~/.codex/config.toml`
- Convex CLI: `1.31.7`
- Date window: 2026-02-18

## Initial symptom

Convex CLI commands worked, but Convex MCP tool calls inside Codex failed with:

```text
Not Authorized: Run `npx convex dev` to login to your Convex project.
```

Observed contradiction:

- `npx convex dev` succeeded and provisioned/connected the dev deployment.
- `npx convex data` succeeded and listed tables.
- `mcp__convex__status` (through Codex MCP) failed with `Not Authorized`.

## MCP config at the time of failure

```toml
[mcp_servers.convex]
command = "npx"
args = [
  "-y", "convex@latest", "mcp", "start",
  "--project-dir", "/home/grace/projects/chem-model-edit",
  "--env-file", "/home/grace/projects/chem-model-edit/.tmp/dev-secrets.env"
]
```

## Working assumptions before deep investigation

1. Maybe Codex had not been restarted after config changes.
2. Maybe `CONVEX_DEPLOY_KEY` in the env file was stale or malformed.
3. Maybe CLI and MCP were authenticating through different paths.
4. Maybe `npx -y convex@latest` introduced instability unrelated to Convex auth itself.

## Investigation log (chronological)

### 1) Verify local secrets and env-file contents safely

We checked key names and non-empty state only (no secret values printed):

- `CONVEX_DEPLOYMENT=SET`
- `CONVEX_URL=SET`
- `CONVEX_SITE_URL=SET`
- `CONVEX_DEPLOY_KEY=SET` (at that stage)

We also confirmed `~/.convex/config.json` existed with `accessToken=SET`.

### 2) Confirm whether auth failure is MCP-only

- `npx convex data --env-file .tmp/dev-secrets.env` succeeded.
- `npx convex data` succeeded.
- `mcp__convex__status` still failed.

Conclusion: failure was specific to MCP path, not a total Convex login failure.

### 3) Inspect Convex MCP/CLI implementation behavior

We reviewed Convex CLI source for MCP/auth flow:

- `mcp start` request handler calls `initializeBigBrainAuth(...)` then `checkAuthorization(...)` per tool call.
- If authorization check returns non-200, MCP tool returns the exact `Not Authorized` message.

Important implementation detail:

- When `--env-file` is provided, Convex reads `CONVEX_DEPLOY_KEY` from that file for auth initialization.
- That key can take precedence over global login token.

### 4) Test removing `CONVEX_DEPLOY_KEY`

We tried removing `CONVEX_DEPLOY_KEY` from `.tmp/dev-secrets.env` while keeping `--env-file` configured.

Result:

- `npx convex env list --env-file .tmp/dev-secrets.env` failed with `MissingAccessToken` / request to authenticate.

Interpretation:

- With `--env-file` path active, missing deploy key can break auth selection for these commands.

### 5) Test placeholder deploy key strategy

We tested a placeholder value for `CONVEX_DEPLOY_KEY` to force fallback.

Result:

- Convex rejected it (`Please set CONVEX_DEPLOY_KEY to a new key...`).

### 6) Identify an additional instability with `-y convex@latest`

In one probe run, `npx -y convex@latest ...` hit npm cache permission issues (`EACCES` in `~/.npm/_cacache`).

Even though this was not the core auth issue, it increased nondeterminism.

### 7) Final mitigation

We changed MCP startup strategy to avoid env-file based auth entirely:

```toml
[mcp_servers.convex]
command = "npx"
args = [
  "convex", "mcp", "start",
  "--project-dir", "/home/grace/projects/chem-model-edit"
]
```

And we kept `.tmp/dev-secrets.env` with only deployment-routing keys:

- `CONVEX_DEPLOYMENT`
- `CONVEX_URL`
- `CONVEX_SITE_URL`

(`CONVEX_DEPLOY_KEY` removed from this file.)

After restarting Codex (so MCP server process restarted), `mcp__convex__status` succeeded.

## Final verification snapshot

`mcp__convex__status` returned available deployment data (own dev deployment), and follow-up calls worked:

- `envList`: succeeded (`variables: []` at that moment)
- `tables`: succeeded (returned table schemas)

## Root-cause statement (practical)

Most likely root cause was an auth-path conflict introduced by `--env-file` + deploy key handling in MCP startup, combined with process-lifecycle confusion (old MCP process surviving config changes).

In plain terms:

- CLI was usable,
- MCP server auth check path still failed,
- and removing env-file-based key resolution plus forcing clean restart made MCP stable.

## Operational guidance

### Recommended baseline for local Codex + Convex MCP

- Prefer:
  - `npx convex mcp start --project-dir <repo>`
- Avoid `-y convex@latest` in MCP config unless you explicitly need floating upgrades.
- Use `--env-file` only when you intentionally need env-file-driven deployment/key selection.
- If MCP auth is flaky, validate whether global token exists at `~/.convex/config.json` and restart Codex fully.

### Safe troubleshooting order

1. Verify CLI: `npx convex data` or `npx convex env list`.
2. Verify MCP: `mcp__convex__status`.
3. If CLI works but MCP fails:
   - inspect MCP startup args,
   - simplify to `npx convex mcp start --project-dir ...`,
   - restart Codex.
4. Reintroduce `--env-file` only if required.

## Security notes

- Do not paste raw `CONVEX_DEPLOY_KEY` in logs/issues/docs.
- Validate key existence and non-empty state without printing value.
- Keep local secrets in `.tmp/dev-secrets.env` as project policy requires.

## References (web + source)

- Convex MCP Server docs (setup/options/tools):
  - https://docs.convex.dev/ai/convex-mcp-server
- Convex CLI docs (`convex dev`, project configuration):
  - https://docs.convex.dev/cli
- Convex deploy key types and behavior:
  - https://docs.convex.dev/cli/deploy-key-types
- Convex CLI MCP implementation (`mcp.ts`):
  - https://unpkg.com/convex@1.31.7/src/cli/mcp.ts
- Convex auth check implementation (`login.ts`):
  - https://unpkg.com/convex@1.31.7/src/cli/lib/login.ts
- Convex deployment/auth initialization (`deploymentSelection.ts`):
  - https://unpkg.com/convex@1.31.7/src/cli/lib/deploymentSelection.ts
- Convex deploy key processing (`utils.ts`):
  - https://unpkg.com/convex@1.31.7/src/cli/lib/utils/utils.ts

## Appendix: commands used during investigation

```bash
# verify CLI path
npx convex data
npx convex env list --env-file .tmp/dev-secrets.env

# inspect MCP option surface
npx convex mcp start --help

# verify CLI version
npx convex --version
```
