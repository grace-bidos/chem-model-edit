# Tasks: One-shot Clerk Migration + Minimal qe.relax.v1

**Input**: `specs/007-clerk-auth-zpe-relax/spec.md` / `specs/007-clerk-auth-zpe-relax/plan.md`  
**Epic**: #151

## Planning PR
- [ ] T700 [P] Link issue hierarchy and finalize plan docs in `specs/007-clerk-auth-zpe-relax/`
- [ ] T701 [P] Add superseded notice and compatibility policy to `specs/004-zpe-user-queue-auth/*`
- [ ] T702 [P] Open draft PR with required links/decisions/review policy

## Implementation PRs (By Child Issue)

### PR-A: #153 API Clerk Auth Foundation
- [ ] T710 [P] Add Clerk JWT JWKS verification dependency to API auth layer
- [ ] T711 [P] Add allowlist enforcement and auth mode settings (`AUTH_MODE`, issuer/audience, allowed emails)
- [ ] T712 [P] Add dev-bypass path for local/CI and tests for auth pass/fail cases

### PR-B: #154 API Contract Cleanup
- [ ] T720 [P] Remove legacy `/api/auth/*` routers/schemas/dependencies
- [ ] T721 [P] Export OpenAPI and regenerate `packages/api-client`
- [ ] T722 [P] Update compile/type usages to remove legacy auth client methods

### PR-C: #155 Web Clerk Migration
- [ ] T730 [P] Integrate Clerk provider/hooks into web app shell
- [ ] T731 [P] Replace ToolPanel legacy sign-in/register flow with Clerk flow
- [ ] T732 [P] Ensure protected API requests attach Clerk Bearer token

### PR-D: #156 Ownership Alignment
- [ ] T740 [P] Migrate job owner identity source to Clerk subject
- [ ] T741 [P] Validate queue target ownership checks under Clerk identity
- [ ] T742 [P] Add owner/non-owner regression tests for status/result/files

### PR-E: #157 Minimal qe.relax.v1
- [ ] T750 [P] Add `qe.relax.v1` request/dispatch handling
- [ ] T751 [P] Add minimal worker-side execution and normalized output
- [ ] T752 [P] Add regression tests ensuring `qe.zpe.v1` remains unaffected

## Optional Research
- [ ] T760 [P] (#158) Produce orchestration comparison memo in `investigations/`

## PR Transparency Rule (Definition of Done Addendum)
- [ ] Every implementation PR must include both `Temporary Behavior` and `Final Behavior` sections.
- [ ] If temporary behavior exists, PR must link a concrete follow-up issue/PR before merge.
- [ ] Child issue cannot be closed until temporary behavior is removed or explicitly accepted in the parent epic.
