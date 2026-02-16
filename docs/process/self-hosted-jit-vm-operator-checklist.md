# Self-hosted JIT+VM Operator Checklist

This checklist is the concrete operator sequence for local rollout.

## 0) One-time GitHub setup (operator action)

1. Create a dedicated runner group for CI only.
2. Restrict the group to this repository.
3. Set workflow variable:
   - `CI_SELF_HOSTED_TRUSTED_ROUTING=false` (start disabled)
4. Confirm required label set:
   - `self-hosted`, `linux`, `x64`, `chem-trusted-pr`

## 1) Auth and token prerequisites (operator action)

Use a principal that can manage self-hosted runners at org or repo scope.

- If using `gh` locally:
  - `gh auth login`
  - verify org/repo admin scope for Actions runner APIs
- If using GitHub App:
  - ensure app can create JIT runner configuration

## 2) VM base image prerequisites (operator action, sudo likely required)

Install required packages on the VM image:

- `curl`
- `jq`
- `tar`
- `ca-certificates`
- `git` (optional but recommended for diagnostics)

Also create a non-root user dedicated to runner execution.

## 3) Request JIT config (controller side)

Example (repo scope):

```bash
scripts/runner/request_jit_config.sh \
  --scope repo \
  --owner grace-bidos \
  --repo chem-model-edit \
  --runner-group-id <RUNNER_GROUP_ID> \
  --labels "self-hosted,linux,x64,chem-trusted-pr" \
  --name-prefix chem-jit \
  --out /tmp/jit-config.json
```

## 4) Launch ephemeral runner on VM

```bash
scripts/runner/run_ephemeral_runner.sh \
  --jit-config /tmp/jit-config.json \
  --runner-home /opt/actions-runner
```

If dependency warning appears on first boot, run once:

```bash
cd /opt/actions-runner
sudo ./bin/installdependencies.sh
```

## 4.1) Repair or bootstrap base runner registration (idempotent)

Use this when a long-lived local runner directory exists but registration was deleted
or became invalid.

Dry-run first:

```bash
scripts/runner/bootstrap_local_runner_registration.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --runner-home /opt/actions-runner \
  --labels "self-hosted,linux,x64,chem-trusted-pr" \
  --dry-run
```

Apply:

```bash
scripts/runner/bootstrap_local_runner_registration.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --runner-home /opt/actions-runner \
  --labels "self-hosted,linux,x64,chem-trusted-pr"
```

Behavior summary:

- checks local `.runner` and GitHub runner inventory by name
- if missing/invalid, fetches a fresh registration token and runs `config.sh --unattended`
- reinstalls and starts `svc.sh` systemd service
- keeps existing work folder and runner name defaults unless explicitly overridden

## 5) Tear-down guarantees

- Enforce `1 VM = 1 Job`.
- Destroy VM immediately after runner exits.
- Do not reuse runner workspace between jobs.

## 6) Enable traffic (canary)

After one successful dry run:

- set `CI_SELF_HOSTED_TRUSTED_ROUTING=true`
- open a trusted PR and verify route job selects self-hosted label target

## 7) Fast rollback

Immediate rollback:

- set `CI_SELF_HOSTED_TRUSTED_ROUTING=false`

Hard rollback:

- remove repository access from the dedicated runner group.

## 8) Local runner health check (incident triage)

Use this to compare local `actions.runner.*` service state against GitHub runner status:

```bash
scripts/runner/check_local_runner_health.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --labels "self-hosted,linux,x64,chem-trusted-pr"
```

Expected result:

- exit code `0`: healthy
- non-zero: degraded (service mismatch and/or offline state)

Optional strict mode:

```bash
scripts/runner/check_local_runner_health.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --strict-gh
```

## 9) One-command local base-runner recovery

When jobs are queued or the local runner is stuck offline, run:

```bash
RUNNER_OWNER=grace-bidos \
RUNNER_REPO=chem-model-edit \
RUNNER_LABELS="self-hosted,linux,x64,chem-trusted-pr" \
RUNNER_GROUP="Default" \
scripts/runner/recover_base_runner.sh --runner-home /opt/actions-runner/actions-runner
```

Dry-run preview (safe):

```bash
RUNNER_OWNER=grace-bidos \
RUNNER_REPO=chem-model-edit \
RUNNER_LABELS="self-hosted,linux,x64,chem-trusted-pr" \
RUNNER_GROUP="Default" \
scripts/runner/recover_base_runner.sh \
  --runner-home /opt/actions-runner/actions-runner \
  --dry-run
```

Required env args are always:

- `RUNNER_OWNER`
- `RUNNER_REPO`
- `RUNNER_LABELS`
- `RUNNER_GROUP`

## 9.1) Scale local runner pool (baseline/max/target)

Use this when you want to run multiple trusted PR jobs in parallel on local self-hosted runners.

Example: keep floor `1`, allow up to `4`, and scale now to `4`.

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --baseline 1 \
  --max 4 \
  --target 4
```

Notes:

- `--target` is optional. If omitted, baseline is applied.
- Effective target is clamped to `[baseline, max]`.
- Runner naming remains `home-self-host`, `home-self-host-2`, ...

## 10) Emergency fallback to hosted routing

If self-hosted runners are unavailable during an incident, immediately route trusted PRs back to hosted:

- set repository variable `CI_SELF_HOSTED_TRUSTED_ROUTING=false`

This is safe to apply before or during local runner recovery.
