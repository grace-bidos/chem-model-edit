# Web App

TanStack Start app (SPA mode) hosted on Cloudflare Workers.

## Local development

```bash
pnpm install
pnpm -C apps/web dev
```

## Runtime contract with API

- Backend API base is `/api` and OpenAPI is served at `/api/openapi.json`.
- `API_BASE_PUBLIC` / `API_BASE` / `VITE_API_BASE` should point to the API origin with `/api` included.
  - Example: `https://chem-model-api-668647845784.asia-northeast1.run.app/api`
- `apps/web/src/server/api.ts` is the only HTTP boundary from web to backend.
- API request/response fields are `snake_case` on the wire.

## Cloudflare Workers settings

- Keep `TSS_SHELL=true` in worker vars (required for SPA shell rendering).
- Set `API_BASE_PUBLIC` in the worker environment to the backend `/api` base.
  - Production example:
    - `API_BASE_PUBLIC=https://chem-model-api-668647845784.asia-northeast1.run.app/api`

## Contract update workflow

Any backend API schema change must update generated client artifacts:

```bash
uv sync --dev --project apps/api
uv run --project apps/api python scripts/export_openapi.py
pnpm -C packages/api-client run generate
git diff -- packages/api-client/openapi/openapi.json packages/api-client/src/generated/schema.ts
```

If `git diff` is non-empty, commit generated artifacts in the same PR.

## Verification commands

```bash
pnpm -C apps/web lint
pnpm -C apps/web typecheck
pnpm -C apps/web test
```
