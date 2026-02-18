# Structure Store Persistence (Convex-First)

This document defines the current persistence contract for parsed/imported structures.

## Goal

- Keep `FastAPI on Modal` stateless.
- Persist application-facing structure data in Convex.
- Keep execution system-of-truth in user-managed AiiDA + Slurm.

## Current Ownership

- Structure import/parse artifacts (`structure`, `cif`, optional raw input) are stored in Convex.
- Runtime execution state remains in Convex projection state + execution backend records.
- Scientific workflow/execution source-of-truth remains AiiDA (PostgreSQL on user-managed environment).
- FastAPI keeps no durable local DB for structure storage.

## Convex Tables

- `structures`
  - Key: `tenant_id + workspace_id + structure_id`
  - Payload: `source`, `cif`, `structure`, `params`, timestamps, chunk count.
- `structure_raw_chunks`
  - Key: `tenant_id + workspace_id + structure_id + chunk_index`
  - Payload: chunk text for large raw input.

## API Service Behavior

`apps/api/services/structures.py`

- `STRUCTURE_STORE_BACKEND=convex|memory`.
- Default is `convex`.
- Convex backend requires:
  - `CONVEX_URL`
  - `CONVEX_DEPLOY_KEY`
- `memory` backend is intended for tests/local isolated runs.
- No SQLite persistence is used for structures.

## Why chunk raw input

Quantum ESPRESSO input can be large. Raw input is chunked and persisted in `structure_raw_chunks` to avoid oversized single-record payloads and to keep retrieval deterministic.

## Viewer Data Loading

Web viewer should not rely on browser direct fetch for protected API routes.

- For `/api/structures/{id}/view?format=cif`, web uses authenticated API helper (`requestApi`) and loads CIF text through the backend contract.
- This avoids cross-origin/auth mismatch from direct unauthenticated fetch paths.

## Migration Notes

- Legacy SQLite-based structure persistence is removed from active path.
- Existing tests run with `STRUCTURE_STORE_BACKEND=memory` fixture by default.
- Production/dev with Convex must explicitly set Convex credentials.
