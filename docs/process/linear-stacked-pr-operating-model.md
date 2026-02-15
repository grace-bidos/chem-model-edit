# Linear + Stacked PR Operating Model

This document defines the default delivery model for this repository.
It standardizes how we plan work in Linear, implement with stacked PRs, and pass CI/review/merge gates.

## 1. Scope and goals

- Reduce review bottlenecks.
- Improve safe parallel delivery.
- Keep planning source-of-truth in Linear.
- Keep GitHub focused on PR review and merge.

## 2. Tool roles

- Linear: planning and execution source-of-truth (projects, cycles, issues, priorities, status).
- Graphite (`gt`): primary stacked PR workflow.
- GitHub (`gh`): PR review, merge, comments, checks visibility.
- Git (`git`): low-level repo operations.
- Jujutsu (`jj`): optional local-history helper only when needed; not used for push/PR submission.

## 3. Planning hierarchy

- Roadmap: optional.
- Project: required for each initiative.
- Cycle: required for active delivery windows.
- Issue: required unit of implementation.
- Fixed implementation hierarchy: `Project -> Parent Capability Issue -> Child Delivery Issue` (one nesting level only).
- Active parent and child issues are in cycle by default.
- Exception: a parent capability issue may stay out of cycle only when intentionally long-running and not active this week.

Current default for this phase:

- No roadmap object.
- Active projects are used as top-level planning containers.
- Cycle duration is 1 week unless explicitly changed by team settings.

## 4. Work item taxonomy

Every implementation issue must have:

- Type axis: `Ask`, `Show`, or `Ship`.
- Size axis: `XS`, `S`, `M`, or `L`.

### 4.1 Type definitions

- Ask: spec/plan/architecture/DB/security-impacting decisions that require agreement first.
- Show: independent feature or behavior changes (including UI or refactor with runtime impact).
- Ship: docs/tests/enablement/polish tasks that finalize readiness.

### 4.2 Size policy

- XS/S: preferred for fast review.
- M: allowed with clear scope boundary.
- L: split before merge queue entry.

## 5. Branch, commit, and PR conventions

- Branch format: `feature/identifier-title`.
- Example: `feature/GRA-21-define-cutover-flags`.
- Commit messages include issue key (example: `feat: add cutover flags (GRA-21)`).
- PR title starts with issue key (example: `GRA-21: Define cutover flags and staged migration`).
- One issue maps to one PR layer in a stack.

## 6. Merge Queue policy

### 6.1 Queue entry rules

- `Ask`:
  - If runtime-impacting (architecture/DB/security/API contract), queue required.
  - If docs-only decision artifact, queue optional.
- `Show`: queue required.
- `Ship`:
  - Docs/tests-only can be queue optional.
  - Runtime-touching changes require queue.

### 6.2 Size gate for queue

- `XS/S/M`: queue eligible.
- `L`: not eligible; split first.

### 6.3 CI gates

Required pre-merge checks (fast):

- lint
- typecheck
- unit/smoke
- policy check (type/size labels or equivalent metadata)

Optional post-merge checks (expand later):

- heavy integration or E2E suites

## 7. End-to-end execution flow

1. Plan in Linear:
   - create/update project and cycle
   - ensure issue has type, size, priority, acceptance criteria
2. Start implementation:
   - create branch from issue key with `gt create` or equivalent
   - keep one issue per branch/PR
3. Build stack:
   - stack related issues in dependency order (Ask -> Show -> Ship)
4. Submit:
   - `gt submit` for current stack scope
5. Review:
   - review in GitHub
   - keep fixes inside same stack branch
   - include architecture/contract context in PR description when runtime boundaries change
6. Merge:
   - merge bottom-up: land the base PR first, then proceed upward through the stack
   - use merge queue where required
7. Sync and continue:
   - `gt sync` and restack as needed
   - advance Linear status (`Todo` -> `In Progress` -> `In Review` -> `Done`)

## 8. Backend refresh mapping (historical example)

This section is a historical snapshot from the initial backend-refresh rollout.
Use it as an example only; current planning source-of-truth is Linear.

### 8.1 Linear projects

- `Backend Refresh (Convex + AiiDA/Postgres + Slurm)`
- `Specialized Services Exploration (Low Priority)`

### 8.2 Epics

- `GRA-10`: backend refresh epic
- `GRA-11`: specialized exploration epic

### 8.3 Cycle loading rule

- Cycle 1 (1 week) included active parent capability `GRA-10` and active backend-refresh child issues (`GRA-12` to `GRA-25`).
- `GRA-37` is a follow-up Ask ADR and is loaded independently when architecture clarification is required.
- Parent `GRA-11` and specialized exploration (`GRA-26` to `GRA-30`) stayed out of Cycle 1 while not active for that week.

### 8.4 Suggested stack slices for backend refresh

- Stack A (Ask/design first):
  - `GRA-12` -> `GRA-13` -> `GRA-17`
- Stack B (core runtime path):
  - `GRA-14` -> `GRA-15` -> `GRA-16` -> `GRA-18` -> `GRA-20`
- Stack C (migration and operations):
  - `GRA-21` -> `GRA-22` -> `GRA-19` -> `GRA-23` -> `GRA-24` -> `GRA-25`
  - follow ADR-0003 single-write cutover policy for migration PRs (`GRA-21`, `GRA-22`): no long-lived dual-write phase, rollback behavior documented, and temporary emergency-stop flags must include explicit expiry/removal plan

### 8.5 Accepted ADRs for backend refresh foundations

- `GRA-12`: `docs/adr/ADR-0001-system-of-record-boundaries.md`
- `GRA-13`: `docs/adr/ADR-0002-jobstate-contract.md`
- `GRA-37`: `docs/adr/ADR-0003-projection-boundary-and-single-write-cutover.md`

### 8.6 Active refactor plan artifacts

- `GRA-14`: `docs/process/aiida-runtime-bootstrap.md`
- `GRA-19`: `docs/process/neon-strategy-validation.md`
- `GRA-20`: `docs/process/fastapi-thin-adapter-cleanup-plan.md`
- `GRA-21`: `docs/process/redis-worker-cutover-flags.md`
- `GRA-22`: `docs/process/redis-worker-retirement-plan.md`
- `GRA-24`: `docs/process/structure-store-persistence.md`

## 9. Review bottleneck controls

- Keep stack depth typically <= 4 PRs when possible.
- Keep each PR focused on one issue only.
- Prefer early draft PR creation to overlap review waiting time with next stacked task.
- Prioritize reviewing lower-stack PRs first to unblock the whole stack.

## 10. Governance

- Linear is source-of-truth for status and priority.
- GitHub Issues are optional mirrors or external-facing discussion threads.
- If mismatch occurs, Linear state wins and GitHub metadata is reconciled.
