# GRA-128: uv.lock dirty delta evaluation (2026-02-17)

## Scope
Compared:
- `/tmp/chem-model-edit-dirty-backup/uv.lock.dirty`
- `apps/api/uv.lock` (current `main`)
- `apps/api/pyproject.toml` (current `main`)

## Findings
- Lockfile drift is very large: `1071 insertions`, `52 deletions`.
- Dirty lock introduces a full `aiida-core` transitive tree (e.g. `aiida-core`, `aio-pika`, `plumpy`, `sqlalchemy`, `psycopg`, etc.).
- Dirty lock also changes versions of already shared packages, including downgrades:
  - `click 8.3.1 -> 8.1.8`
  - `psutil 7.1.3 -> 6.1.0`
  - `tzdata 2025.2 -> 2021.3`
- `apps/api/pyproject.toml` already contains an `aiida` dependency group (`aiida-core>=2.7.3`), but there is no paired backup change indicating deliberate dependency-spec edits for this lockfile expansion.

## Decision
**Discard the dirty lockfile delta** (do not adopt into `apps/api/uv.lock`).

## Rationale
- The delta is not minimal and affects shared dependency versions outside a narrowly scoped, intentional change.
- No explicit dependency declaration change accompanies this lock expansion in the provided backup artifacts.
- Adopting this as-is would introduce high regression risk for API/runtime environments without clear product intent.

## Follow-up (if AIIDA group locking is required)
- Create a focused follow-up change that:
  - documents intended `uv` group-lock policy,
  - updates lock generation command/process explicitly,
  - validates with targeted API test/type/lint checks before merge.
