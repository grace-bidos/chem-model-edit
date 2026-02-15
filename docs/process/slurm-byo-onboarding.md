# Slurm BYO Compute Onboarding (GRA-15 Slice)

This document defines the first practical onboarding/validation slice for integrating Slurm-backed BYO compute nodes with an explicit fallback policy.

Scope of this slice (`GRA-15`):

- Define onboarding requirements for Slurm-backed node registration.
- Define queue/account/QoS mapping policy.
- Provide an executable validation helper plus example config files.
- Keep runtime behavior unchanged (operator prep and validation only).

Related baseline:

- `docs/process/aiida-runtime-bootstrap.md` (AiiDA runtime bootstrap safety model)

## Onboarding Requirements

Before registering a new BYO node, operators must provide:

1. Slurm cluster identity
   - `cluster_name` that matches the registration payload.
2. Scheduler allow-lists
   - Allowed `partitions`, `accounts`, and `qos` values.
3. Logical queue mappings
   - Product-facing queue names mapped to concrete Slurm `partition/account/qos`.
4. Fallback policy
   - Behavior when an unknown queue is requested (`deny` or `route-default`).
5. Registration requirements
   - Required node labels and required scheduler health checks.

Use these artifacts as the starting point:

- Policy config template: `apps/api/config/slurm/onboarding.policy.example.json`
- Registration template: `apps/api/config/slurm/node-registration.example.json`
- Validator script: `scripts/validate-slurm-onboarding.py`

## Queue/Account/QoS Mapping Policy

Each logical queue entry in `queue_mappings` must:

- Use a unique `queue` name.
- Reference only allow-listed `partition`, `account`, and `qos`.
- Set a positive `max_walltime_minutes`.

Fallback policy options:

1. `deny`
   - Unknown queue requests fail validation.
   - Recommended when strict tenant isolation is required.
2. `route-default`
   - Unknown queue requests are redirected to `fallback_policy.default_queue`.
   - `default_queue` must exist in `queue_mappings`.
   - Recommended for phased onboarding where not all clients are upgraded yet.

## Validation Checklist for New Node Registration

Run through this checklist before activation:

1. Policy file validates with zero errors.
2. Cluster name in registration matches policy cluster name.
3. `enabled_queues` only contain mapped queues.
4. All required node labels are present and non-empty.
5. All required health checks are present and `true`.
6. If `deny_if_drain_or_down=true`, node state is not `down`, `drain`, `draining`, or `fail`.
7. Optional queue-resolution smoke test matches expected fallback behavior.

## Validation Commands

Validate policy only:

```bash
python3 scripts/validate-slurm-onboarding.py \
  --policy apps/api/config/slurm/onboarding.policy.example.json
```

Validate policy + node registration + fallback behavior:

```bash
python3 scripts/validate-slurm-onboarding.py \
  --policy apps/api/config/slurm/onboarding.policy.example.json \
  --registration apps/api/config/slurm/node-registration.example.json \
  --requested-queue experimental
```

Success returns exit code `0`; validation errors return exit code `1`.
