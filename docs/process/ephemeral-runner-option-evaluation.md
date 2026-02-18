# Ephemeral Runner Option Evaluation

This document closes the remaining delivery criteria for runner-option validation.

Related issues:

- GitHub: #362
- Linear: GRA-99

## Candidate options

### Option A: ARC (Actions Runner Controller) on Kubernetes

- Summary: Use GitHub ARC scale sets to manage ephemeral runners on K8s.
- Best fit: teams already operating Kubernetes clusters with SRE support.

### Option B: Cloud-managed ephemeral runners

- Summary: Use hosted ephemeral runner products/services from cloud providers.
- Best fit: teams optimizing for low operations burden over maximum control.

### Option C: Local VM pool + JIT registration (current path)

- Summary: Keep current self-hosted controller and JIT workflow with VM pool.
- Best fit: teams needing strict local control and predictable cost envelopes.

## Comparison matrix

| Criterion | Option A ARC/K8s | Option B Cloud-managed | Option C Local VM+JIT |
| --- | --- | --- | --- |
| Setup complexity | High (K8s, ARC lifecycle, RBAC) | Low to medium | Medium |
| Security posture | Strong when cluster hardening is mature | Good defaults, less control | Strong local control; more operator responsibility |
| Cost model | Cluster and ops overhead; efficient at scale | Potentially highest unit cost | Lowest infra unit cost, higher operator time |
| Scaling latency | Good with warm nodes; variable when cluster scales | Good to excellent | Good with warm pool; slower from cold boot |
| Maintenance burden | High | Low | Medium |

## Recommendation

Recommend **Option C (Local VM + JIT)** as the current production path.

Why:

- Already implemented and validated in this repository's CI and operator playbooks.
- Provides direct control over routing and incident recovery.
- Lowest migration risk for current stack and team capacity.

## Rollback / fallback strategy

Primary rollback:

- Disable trusted self-hosted routing immediately:
  - `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body false`

Operational fallback:

- Keep GitHub-hosted fallback lane always available.
- Re-enable trusted routing only after runner health is stable.

Strategic fallback:

- If maintenance burden increases, phase toward Option A (ARC) after formal capacity review.

## Decision record

- Decision date: 2026-02-18
- Decision owner: Repository maintainers
- Revisit trigger:
  - sustained runner incident rate above defined SLO in operator policy, or
  - materially increased trusted PR throughput requirement beyond local VM capacity
