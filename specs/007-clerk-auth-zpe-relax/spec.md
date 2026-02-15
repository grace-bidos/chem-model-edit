# Feature Specification: One-shot Clerk Migration + Minimal qe.relax.v1

**Feature Branch**: `[stack/clerk-plan-pr]`  
**Created**: 2026-02-11  
**Status**: Draft  
**Input**: Epic #151 and child issue #152

## Context

Current auth/session is implemented as local account registration and Redis-backed sessions (`/api/auth/*`).
For this project phase, the operation remains low-traffic and mostly single-user, with occasional testers.
The team decided to replace legacy auth in one shot and standardize on Clerk JWT verification.

## Decisions (Locked)

- Auth migration strategy: **one-shot replacement** (no long-term dual operation)
- API auth mode: **Clerk JWT verification via JWKS**
- Signup policy: **allowlist-only** for low-user research operation
- Local/CI mode: **dev-bypass** is allowed; production remains Clerk-only
- Compute scope for this epic: **`qe.zpe.v1` + minimal `qe.relax.v1`**
- Orchestrator future options (Prefect/Temporal/FireWorks/AiiDA): **research issue only**, no implementation in this epic

## User Scenarios & Tests (Required)

### User Story 1 - Clerk-authenticated users can access protected APIs (Priority: P1)

As an authenticated user, I can call protected API endpoints with Clerk Bearer JWT.

**Independent Test**: A valid Clerk JWT can access `/api/zpe/*` endpoints that require auth.

**Acceptance Scenarios**:

1. **Given** a valid Clerk JWT, **When** `/api/zpe/targets` is called, **Then** 200 is returned.
2. **Given** an invalid/expired JWT, **When** the same endpoint is called, **Then** 401 is returned.
3. **Given** a valid JWT with a non-allowlisted email, **When** the endpoint is called, **Then** 403 is returned.

---

### User Story 2 - Legacy `/api/auth/*` contract is removed (Priority: P1)

As a maintainer, I want legacy auth endpoints removed so the API contract has a single auth model.

**Independent Test**: OpenAPI and generated API client no longer include `/api/auth/register|login|logout|me`.

**Acceptance Scenarios**:

1. **Given** regenerated OpenAPI, **When** client types are generated, **Then** legacy auth models are absent.
2. **Given** web build/typecheck, **When** auth calls are compiled, **Then** no references to legacy auth client methods remain.

---

### User Story 3 - Web auth flow is Clerk-native (Priority: P1)

As a user, I use Clerk sign-in and the web app sends Clerk tokens on API requests.

**Independent Test**: Signed-in state from Clerk enables ZPE job submission and target management.

**Acceptance Scenarios**:

1. **Given** a signed-in Clerk user, **When** `createEnrollToken` or `createZpeJob` is called, **Then** request includes Clerk Bearer token.
2. **Given** signed-out state, **When** protected calls are made, **Then** UI shows auth-required state.

---

### User Story 4 - ZPE ownership checks remain correct with Clerk subject (Priority: P1)

As a user, I can only see my own job status/results after migration.

**Independent Test**: Job owner access control still blocks non-owners with 403 after identity source changes to Clerk subject.

**Acceptance Scenarios**:

1. **Given** owner submits job, **When** owner fetches status/result/files, **Then** data is returned.
2. **Given** different authenticated user, **When** fetching same resources, **Then** 403 is returned.

---

### User Story 5 - Minimal qe.relax.v1 works on existing worker path (Priority: P2)

As a maintainer, I can run a minimal relax job using the existing queue/worker architecture.

**Independent Test**: `qe.relax.v1` can enqueue, execute, and return normalized result shape without breaking ZPE flows.

**Acceptance Scenarios**:

1. **Given** valid relax payload, **When** job is submitted, **Then** worker executes and result is retrievable.
2. **Given** existing ZPE flow, **When** regression tests run, **Then** ZPE behavior remains unchanged.

## Requirements (Required)

### Functional Requirements

- **FR-001**: API must validate Clerk JWT via JWKS for protected endpoints.
- **FR-002**: API must enforce allowlist-based signup/usage policy.
- **FR-003**: API must provide `AUTH_MODE=clerk|dev-bypass`; production uses `clerk`.
- **FR-004**: Legacy `/api/auth/register`, `/api/auth/login`, `/api/auth/logout`, `/api/auth/me` must be removed.
- **FR-005**: Generated API client must reflect the updated OpenAPI contract with no legacy auth methods.
- **FR-006**: Web app must use Clerk provider/hooks; legacy local session storage flow is removed from primary path.
- **FR-007**: Job ownership and queue-target authorization must map to Clerk subject identity.
- **FR-008**: Existing worker token flow (`/api/zpe/compute/*`) must remain compatible.
- **FR-009**: Minimal `qe.relax.v1` must be supported on the existing queue/worker path.

### Non-Functional Requirements

- **NFR-001**: Keep Redis-based queue/worker architecture in this epic.
- **NFR-002**: Keep PR size small and child-issue scoped to reduce CodeRabbit wait/rate-limit impact.
- **NFR-003**: Design comments from CodeRabbit must be resolved in both spec and plan docs.

## Public API / Interface Changes

- Remove legacy auth endpoints:
  - `POST /api/auth/register`
  - `POST /api/auth/login`
  - `POST /api/auth/logout`
  - `GET /api/auth/me`
- Protected APIs continue to use `Authorization: Bearer <token>`, now with Clerk JWT.
- Optional diagnostics endpoint can be added as `GET /api/me`.
- Backend env contract:
  - `AUTH_MODE=clerk|dev-bypass`
  - `AUTH_CLERK_ISSUER`
  - `AUTH_CLERK_AUDIENCE` (optional)
  - `AUTH_ALLOWED_EMAILS`

## Success Criteria

- **SC-001**: Planning PR is merged with no unresolved high-impact design feedback from CodeRabbit.
- **SC-002**: Child issues have decision-complete acceptance criteria and fixed implementation order.
- **SC-003**: Implementation can proceed with no additional architecture-level decisions.

## Out Of Scope

- Prefect/Temporal/FireWorks/AiiDA integration implementation
- CP2K support implementation
- Queue backend replacement (Redis -> other)
