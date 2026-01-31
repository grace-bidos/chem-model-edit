# Epic: Editor UI Replace (editor)

## Goal
Build a new editor experience at `/editor` and iterate toward the final UI. The new UI should center on **immutable structures**: imported `.in` data is treated as read-only, and tools generate new structures instead of in-place editing.

## Principles
- **Immutable inputs**: imported files are preserved; no direct edits to original data.
- **Generated outputs**: tools (Transfer / Supercell / etc.) create new structures.
- **Provenance-first**: each structure should record how it was created (source + operation).
- **Low-risk rollout**: legacy editor can remain unchanged; parity is not required.

## Scope
- New UI layout and navigation at `/editor`.
- File Manager lists structures + creation lineage.
- Workspace panels render per-structure views.
- Tool panels provide structure-generating actions.
- ZPE/Vibrations stays **preview-only** for now.

## Non-goals (for this epic)
- Full parity with the legacy editor.
- Rewriting backend APIs.
- ZPE/Vibrations feature wiring.

## Definition of Done
- `/editor` is usable and matches the new mental model.
- Structure creation flows are centered on **generation** instead of direct edits.
- The legacy editor does not need full functional parity.

## Progress Log
- 2026-01-10: `/editor` route + mock layout + tool renaming + vibration preview shell.

## Notes
Feature Issues will be created only as needed during UI review and iteration.
