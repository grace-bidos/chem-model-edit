# Codex Collaboration Sub-Agent Pilot (Backend Refresh)

This document defines a 1-day pilot where one primary Codex session orchestrates focused sub-agents to accelerate Backend Refresh work without losing safety or ownership clarity.

## 1. Goals

- Increase same-day throughput for Backend Refresh tasks.
- Reduce review wait by parallelizing independent implementation lanes.
- Keep decisions, status, and acceptance criteria aligned with Linear and stacked PR flow.
- Preserve quality gates (`lint`, `typecheck`, `unit/smoke`) while moving faster.

## 2. Lane Structure (4+1)

Use five concurrent lanes: four execution lanes plus one control lane.

- Lane 1: API contracts and router/service updates.
- Lane 2: Worker/job orchestration and queue behavior.
- Lane 3: Data and migration concerns (schemas, compatibility, safety checks).
- Lane 4: Web/API client integration and contract-consumer updates.
- Lane +1 (Control, Primary Codex): orchestration, dependency management, integration review, and final merge readiness.

## 3. Ownership Rules

- Primary Codex owns scope slicing, lane assignment, integration plan, and final acceptance decision.
- Each sub-agent owns only its assigned lane deliverables and must not change lane boundaries without handoff approval.
- One implementation issue maps to one PR layer; cross-lane dependencies are stacked explicitly.
- Primary Codex is the only actor that resolves cross-lane conflicts and updates the final status in Linear.

## 4. Sub-Agent Usage Boundaries

- Sub-agents may implement code, tests, and local docs inside their lane (for example README updates or short design notes tied to the lane change).
- Sub-agents must not perform cross-lane refactors, policy changes, or unapproved architecture changes.
- Sub-agents must not merge PRs or rewrite other lanes’ commit history.
- Any security/auth, API contract break, or migration risk escalation is handed back to primary Codex immediately.

## 5. Daily Cadence (1-Day Pilot)

1. 09:00-09:20: Primary Codex defines lane briefs (goal, files, acceptance checks, dependency notes).
2. 09:20-12:00: Sub-agents execute first pass; primary Codex runs asynchronous blocker triage.
3. 12:00-12:20: Midday sync; confirm scope, reassign blocked items, freeze risky expansion.
4. 12:20-16:00: Second pass implementation and tests; prepare stacked PR layers.
5. 16:00-17:00: Primary Codex integration review, final checks, pilot outcome summary.

## 6. Risk Controls

- Scope lock: no new lane creation during the pilot day.
- Interface lock: cross-lane contract changes require explicit primary approval.
- Check gate: do not mark lane complete with failing required checks.
- Conflict gate: all merge/rebase conflict resolution is done by primary Codex.
- Rollback readiness: keep each lane deliverable independently revertible through isolated commits and PR layers.

## 7. Pilot Checklist

- Pilot issues in Linear are tagged and mapped to lanes.
- Each lane has written acceptance criteria and file ownership hints.
- Dependency order is defined for stacked PR submission.
- The primary Codex confirms dependency order at pilot start and re-confirms it during the midday sync.
- Required checks are confirmed and runnable locally or in CI.
- End-of-day integration owner (primary Codex) is explicitly assigned.
- Post-pilot retro slot (15-20 min) is booked the same day.

## 8. Success Metrics (1-Day)

- Throughput: at least 4 lane-scoped deliverables reach review-ready state.
- Lead time: median time from lane start to review-ready is under 6 hours.
- Quality: zero known contract regressions and zero merge-with-red-check events.
- Coordination: all blocker escalations acknowledged by primary Codex within 30 minutes.
- Predictability: at least 80% of planned pilot scope is completed or cleanly restacked.

## 9. Pilot Exit Decision

At day end, keep this model for Backend Refresh only if all are true:

- Success metrics are met or only one metric misses by a small margin.
- No severe integration incident occurred.
- Team retro confirms the +1 control lane reduced confusion rather than adding overhead.

## 10. Small-Scale Validation Run (This Repository)

This pilot model was validated with a low-risk trial before full rollout.

- Scope: create one new process document under `docs/process/` only.
- Ownership: one sub-agent owned exactly one file; primary Codex reviewed and accepted output.
- Boundary check: no existing files were modified by the sub-agent.
- Quality check: markdown lint command was not available in this environment, so lint was skipped and recorded.
- Result: trial succeeded as a safe rehearsal for lane-based delegation and primary-controlled integration.
- Runtime note: updating `AGENTS.md` during an active session did not retroactively change sub-agent output language in this trial; treat language-rule updates as effective from a fresh session and keep explicit spawn-level language instructions as a fallback.
