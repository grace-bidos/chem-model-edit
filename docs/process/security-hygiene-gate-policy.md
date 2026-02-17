# Security and Hygiene Gate Policy (Phase-1 -> Phase-2)

This document defines how we run security/hygiene checks while keeping CI phase-1 non-blocking, and how we promote those checks to required gates.

## Scope

- API lane checks in CI:
  - `bandit`
  - `pip-audit`
  - `gitleaks`
  - `deptry`
  - `vulture`
  - `mutmut` smoke
- Local developer entrypoints in `Justfile`:
  - `api-security`
  - `api-security-bandit`
  - `api-security-audit`
  - `api-deadcode`
  - `api-quality-security-hygiene`
  - `api-quality-phase1`

## Phase-1 Operating Mode (Non-blocking)

- CI behavior:
  - Security/hygiene checks run on relevant API changes.
  - Checks are configured as non-blocking (`continue-on-error: true`) in phase-1.
  - Bandit and pip-audit always emit machine-readable reports and a run summary in the job step summary.
  - Reports are uploaded as artifacts for traceability.
- Reviewer expectation:
  - A non-blocking failure is still triaged in PR review.
  - Any high-severity finding is treated as stop-and-fix unless explicitly waived with owner sign-off.

## Promotion Criteria to Phase-2 (Required Gate)

Promotion from phase-1 to required checks is done per tool and must satisfy all criteria below:

- Data window:
  - Observe at least 10 runs on `main` where the tool executed.
- Stability threshold:
  - Last 5 consecutive `main` runs meet baseline with no regressions.
- Baseline threshold:
  - Bandit: `SEVERITY.HIGH == 0`.
  - pip-audit: vulnerability count `== 0`.
  - gitleaks: no active secrets detected.
  - deptry/vulture: no newly introduced findings in touched modules.
- Operational readiness:
  - Findings are reproducible locally via documented `just` commands.
  - Remediation owner is identified for any remaining debt.

If any criterion fails, remain in phase-1 and keep collecting data.

## CI Signal Contract

- CI must publish enough context to make promotion measurable:
  - counts for Bandit and pip-audit in job summary
  - downloadable artifacts for raw reports
  - promotion evidence snapshot in markdown + machine-readable JSON (`.ci-security-promotion-report.md` / `.ci-security-promotion-report.json`)
- Promotion decisions should be recorded in the relevant Linear child issue with run evidence links.

## Promotion Evidence Command

- Local/CI report generation command:
  - `just api-security-promotion-evidence`
- Script:
  - `scripts/gh/security_phase1_promotion_report.py`
- Data contract:
  - reads workflow run/job data via `gh` JSON APIs
  - evaluates promotion criteria over a 10-run measured window
  - checks latest 5 consecutive runs for baseline stability
  - emits markdown summary and JSON payload for automation
- Exit semantics:
  - default mode is non-blocking for findings/policy misses (exit `0`)
  - non-zero exits are reserved for script/runtime failures
  - use `--strict` to make unmet promotion criteria exit non-zero

## Local Workflow

Recommended local sequence before PR update:

```bash
just api-security-bandit
just api-security-audit
just api-quality-security-hygiene
just api-security-promotion-evidence
```

For broader phase-1 parity checks:

```bash
just api-quality-phase1
```

## Non-blocking Caveat

Phase-1 non-blocking means merge is not mechanically blocked by these checks. It does not mean findings are optional. Teams should treat severe findings as release risk and fix or explicitly waive before merge.
