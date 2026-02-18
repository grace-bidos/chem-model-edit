# API Contract Operations

This project separates frontend and backend responsibilities at the OpenAPI contract boundary.

## Fixed contract rules

- API base path: `/api`
- OpenAPI endpoint: `/api/openapi.json`
- JSON field naming: `snake_case`
- Web client source of truth: generated artifacts in `packages/api-client`

## Environment variables

### FastAPI on Modal

Set CORS allow list and preview regex in the Modal-deployed FastAPI service:

- `CORS_ALLOW_ORIGINS`
  - Comma-separated exact origins.
  - Example:
    - `http://localhost:3000,http://localhost:3001,https://chem-model-edit.tadashi240312.workers.dev`
- `CORS_ALLOW_ORIGIN_REGEX`
  - Regex for preview deployments.
  - Example:
    - `^https://[a-z0-9-]+-chem-model-edit\\.tadashi240312\\.workers\\.dev$`

### Cloudflare Workers (Web)

Set the backend base including `/api`:

- `API_BASE_PUBLIC=https://<workspace>--chem-model-edit-api-api.modal.run/api`
- `MODAL_PROXY_KEY=<proxy key>` (required when API uses Modal proxy auth)
- `MODAL_PROXY_SECRET=<proxy secret>` (required when API uses Modal proxy auth)

`TSS_SHELL=true` must remain enabled to render the SPA shell correctly.

## Required contract sync steps

When backend schemas/routes change:

```bash
uv sync --dev --project apps/api
PYTHONPATH=apps/api uv run --project apps/api python apps/api/scripts/export_openapi.py
pnpm -C packages/api-client run generate
git diff --exit-code packages/api-client/
```

CI enforces these steps in the `contract` job.
