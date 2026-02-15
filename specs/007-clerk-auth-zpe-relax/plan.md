# Implementation Plan: One-shot Clerk Migration + Minimal qe.relax.v1

**Branch**: `[stack/clerk-plan-pr]` | **Date**: 2026-02-11 | **Spec**: `specs/007-clerk-auth-zpe-relax/spec.md`  
**Input**: Epic #151 + child issues #152-#157

## Summary

This planning PR defines the issue hierarchy and implementation sequence for one-shot Clerk migration.
The migration removes legacy local auth endpoints, standardizes JWT verification via Clerk JWKS,
and keeps the current Redis-based worker architecture while adding minimal `qe.relax.v1`.

## Key Decisions

- One-shot migration to Clerk (legacy auth removed, no dual-run in production)
- Clerk JWT verified by API through JWKS
- Allowlist-based user policy
- dev-bypass only for local/CI
- Implementation delivered through child-issue-scoped small PRs

## Technical Context

- Frontend: TanStack Start (`apps/web`)
- Backend: FastAPI (`apps/api`)
- Queue/worker: Redis + existing ZPE remote worker flow
- API contract package: `packages/api-client`
- CI/Review: GitHub Actions + CodeRabbit auto review

## Project Structure

```plaintext
specs/007-clerk-auth-zpe-relax/
├── spec.md
├── plan.md
└── tasks.md
```

## Work Breakdown (Child-Issue Driven)

1. **#152 Plan Docs**
   - Add new `specs/007-*` docs
   - Add superseded note to `specs/004-*`
2. **#153 API Clerk Auth Foundation**
   - Add Clerk JWT verification dependency (JWKS)
   - Add allowlist policy and `AUTH_MODE`
3. **#154 API Contract Cleanup**
   - Remove legacy `/api/auth/*`
   - Regenerate OpenAPI/client
4. **#155 Web Clerk Migration**
   - Replace legacy login/register UI/session flow
   - Inject Clerk token for protected API requests
5. **#156 ZPE Ownership Alignment**
   - Use Clerk subject for owner/queue checks
   - Keep worker token flow unchanged
6. **#157 Minimal `qe.relax.v1`**
   - Add minimal payload handling and worker execution path
   - Add regression protections for ZPE

## PR Strategy

- Parent issue: #151
- Child issues: #152-#157
- Optional research issue: #158
- Merge strategy: merge commit
- PR size target: small, child-issue scoped

## CodeRabbit Feedback Handling Policy

- Categorize comments into:
  - `design-change` (architecture or contract impact)
  - `wording/docs` (text-level quality)
- For `design-change`, update both:
  - `specs/007-clerk-auth-zpe-relax/spec.md`
  - `specs/007-clerk-auth-zpe-relax/plan.md`
- Post an English PR comment summarizing changes made per addressed review thread.
- Exit condition for planning PR:
  - No unresolved high-impact review comments
  - No undecided architecture-level choices

## Implementation Gate

Implementation starts only when all are true:

- Planning PR merged
- Parent issue links to final merged plan
- Child issue acceptance criteria are fixed
- Child issue implementation order is fixed

## Risks & Mitigations

- **Risk**: Auth contract churn impacts web/api-client simultaneously.
  - **Mitigation**: Sequence #153 -> #154 -> #155 and validate typecheck/test at each step.
- **Risk**: Ownership checks regress during identity migration.
  - **Mitigation**: Add explicit owner/non-owner test coverage in #156.
- **Risk**: Scope creep from orchestration discussion.
  - **Mitigation**: Track only in #158, no implementation in this epic.
