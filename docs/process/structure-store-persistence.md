# Structure Store Persistence Guide (`GRA-24`)

This document describes the local durable storage model for structure records
and the reset/backup procedure for development.

## Overview

`services.structures` now supports a SQLite-backed store so imported structure
records survive API process restarts.

- backend selector: `STRUCTURE_STORE_BACKEND`
- supported backends: `sqlite` (default), `memory`
- SQLite path: `STRUCTURE_STORE_DB_PATH`
- reset switch: `STRUCTURE_STORE_RESET_ON_START`

Default SQLite path is `<repo>/.just-runtime/structures.sqlite3`.

## Development Defaults

Recommended `.env` values for local restart-survival behavior:

```env
STRUCTURE_STORE_BACKEND=sqlite
STRUCTURE_STORE_DB_PATH=.just-runtime/structures.sqlite3
STRUCTURE_STORE_RESET_ON_START=0
```

## Intentional Reset Procedure (Dev)

Use one of the following options.

Option A (one-time reset on next start):

```bash
STRUCTURE_STORE_RESET_ON_START=1 uv run uvicorn apps.api.main:app --reload --port 8000
```

Option B (manual file delete):

```bash
rm -f .just-runtime/structures.sqlite3
```

## Backup Procedure (Dev/Staging)

SQLite copy backup:

```bash
cp .just-runtime/structures.sqlite3 .just-runtime/structures.sqlite3.bak
```

Restore:

```bash
cp .just-runtime/structures.sqlite3.bak .just-runtime/structures.sqlite3
```

## Notes

- ZPE queue/lease/result paths remain separate concerns and are not changed by
  this structure-store persistence layer.
- Tests for restart survival live in
  `apps/api/tests/test_structures_persistence.py`.
