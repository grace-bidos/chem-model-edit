# Web App

TanStack Start app (SPA mode) hosted on Cloudflare Workers.

## Local development

```bash
pnpm install
pnpm -C apps/web dev
```

## Runtime contract with API

- Backend API base is `/api` and OpenAPI is served at `/api/openapi.json`.
- Runtime variable precedence:
  - `API_BASE_PUBLIC`: preferred in Cloudflare Workers runtime (server-side).
  - `API_BASE`: server-side fallback (Node/SSR environments).
  - `VITE_API_BASE`: build-time fallback for local/dev.
- All of the above should point to the API origin with `/api` included.
  - Example: `https://<workspace>--chem-model-edit-api-api.modal.run/api`
- `apps/web/src/server/api.ts` is the only HTTP boundary from web to backend.
- API request/response fields are `snake_case` on the wire.

## Cloudflare Workers settings

- Keep `TSS_SHELL=true` in worker vars (required for SPA shell rendering).
- Set `API_BASE_PUBLIC` in the worker environment to the backend `/api` base.
  - Production example:
    - `API_BASE_PUBLIC=https://<workspace>--chem-model-edit-api-api.modal.run/api`

## Contract update workflow

Any backend API schema change must update generated client artifacts:

```bash
uv sync --dev --project apps/api
PYTHONPATH=apps/api uv run --project apps/api python apps/api/scripts/export_openapi.py
pnpm -C packages/api-client run generate
git diff -- packages/api-client/
```

If `git diff` is non-empty, commit generated artifacts in the same PR.

## SPA Shell (Production)

This app uses TanStack Start in SPA mode. In production, the worker must run in
shell mode so the root HTML + scripts are rendered. Ensure `TSS_SHELL=true`
is set in the deploy environment (already configured in `wrangler.jsonc` for
Cloudflare Workers).

## Verification commands

```bash
pnpm -C apps/web lint
pnpm -C apps/web typecheck
pnpm -C apps/web test
pnpm -C apps/web test:coverage
pnpm -C apps/web test:a11y
pnpm -C apps/web test:fastcheck
pnpm -C apps/web test:mutation
pnpm -C apps/web test:e2e:install
pnpm -C apps/web test:e2e
pnpm -C apps/web test:e2e:smoke
```

## Local-first quality workflow

- CI required gates are intentionally minimal: `lint`, `typecheck`, `test`, `a11y`.
- Heavier checks are local-first:
  - Mutation testing: `pnpm -C apps/web test:mutation`
  - Visual verification: `pnpm -C apps/web build-storybook` and `pnpm -C apps/web chromatic`
- Optional E2E smoke: `pnpm -C apps/web test:e2e:smoke`
- Playwright is local-only for now (not CI required). Run `pnpm -C apps/web test:e2e:install` before first run.
- Strictness and promotion rubric: `docs/process/web-quality-playbook-v2.md`.
