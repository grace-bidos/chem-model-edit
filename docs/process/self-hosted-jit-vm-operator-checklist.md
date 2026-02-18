# Self-hosted JIT+VM Operator Checklist

This checklist is the concrete operator sequence for local rollout.

## Standard policy (quick reference)

- Supervisor mode is primary for steady operations:
  - `scripts/runner/setup_pool_supervisor_one_command.sh`
- Trusted routing default is enabled:
  - `CI_SELF_HOSTED_TRUSTED_ROUTING=true`
- Timer reconcile mode is fallback only:
  - `scripts/runner/setup_pool_reconcile_one_command.sh`
- Fast helper wrappers:
  - token refresh setup: `scripts/runner/setup_github_app_token_refresh_one_command.sh`
  - base-runner recovery: `scripts/runner/recover_base_runner_one_command.sh`
  - guarded routing rollback: `scripts/runner/guard_trusted_routing.sh`

## 0) One-time GitHub setup (operator action)

1. Create a dedicated runner group for CI only.
2. Restrict the group to this repository.
3. Set workflow variable:
   - `CI_SELF_HOSTED_TRUSTED_ROUTING=false` during initial setup
4. Confirm required label set:
   - `self-hosted`, `linux`, `x64`, `chem-trusted-pr`

## 1) Auth and token prerequisites (operator action)

Use a principal that can manage self-hosted runners at org or repo scope.

- If using `gh` locally:
  - `gh auth login`
  - verify org/repo admin scope for Actions runner APIs
- If using GitHub App:
  - ensure app can create JIT runner configuration
  - prefer installation access tokens (short-lived) for runner operations
  - keep app private key outside repository and host it with restricted file permissions

Recommended default:

- Use GitHub App installation tokens for routine runner operations.
- Keep PAT flow only for emergency break-glass recovery.

Helper commands for GitHub App setup:

```bash
scripts/runner/check_github_app_jwt.sh \
  --app-id <APP_ID> \
  --private-key-file <PRIVATE_KEY_PATH>

scripts/runner/get_github_app_installation_id.sh \
  --app-id <APP_ID> \
  --private-key-file <PRIVATE_KEY_PATH> \
  --owner grace-bidos \
  --repo chem-model-edit
```

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
  --token-source env \
  --token-env-var GH_TOKEN \
  --out /tmp/jit-config.json
```

When using GitHub App tokens, export `GH_TOKEN` with a short-lived installation token first:

```bash
export GH_TOKEN="$(
  scripts/runner/request_github_app_token.sh \
    --app-id <APP_ID> \
    --installation-id <INSTALLATION_ID> \
    --private-key-file <PRIVATE_KEY_PATH>
)"
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

## 6) Validate traffic routing (canary)

After setup + one successful dry run, enable trusted routing and verify behavior:

- `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body true`
- open a trusted PR and verify route job selects self-hosted label target

Canary policy:

- Stage 0 (smoke): 1 trusted PR
- Stage 1 (light): next 10 trusted PRs
- Stage 2 (broader): next 25 trusted PRs
- Stage 3 (steady state): all trusted PRs

Advance only if rollback triggers are not hit.

## 7) Fast rollback

Immediate rollback:

- `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body false`
- or run guarded rollback helper:
  - `scripts/runner/guard_trusted_routing.sh --owner <owner> --repo <repo>`

Hard rollback:

- remove repository access from the dedicated runner group.

Restore default routing after incident:

- `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body true`

Rollback triggers (execute immediate rollback if any is true):

- two consecutive trusted PRs fail before test execution because runner boot/registration failed
- p95 trusted queue wait exceeds 120 seconds across latest 10 trusted PRs
- self-hosted runner offline state persists for more than 10 minutes after recovery command
- route mismatch is observed (trusted PR not routed to self-hosted labels, or untrusted routed to self-hosted)

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
scripts/runner/recover_base_runner_one_command.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --labels "self-hosted,linux,x64,chem-trusted-pr" \
  --group "Default" \
  --runner-home /opt/actions-runner
```

Dry-run preview (safe):

```bash
scripts/runner/recover_base_runner_one_command.sh \
  --owner grace-bidos \
  --repo chem-model-edit \
  --labels "self-hosted,linux,x64,chem-trusted-pr" \
  --group "Default" \
  --runner-home /opt/actions-runner \
  --dry-run
```

Manual fallback (`recover_base_runner.sh`) required env args:

- `RUNNER_OWNER`
- `RUNNER_REPO`
- `RUNNER_LABELS`
- `RUNNER_GROUP`

## 9.1) Scale local runner pool (min/max/target)

Use this when you want to run multiple trusted PR jobs in parallel on local self-hosted runners.

Example: keep floor `1`, allow up to `4`, and scale now to `4`.

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --min 1 \
  --max 4 \
  --target 4
```

Notes:

- `--target` is optional. If omitted, min is applied.
- Effective target is clamped to `[min, max]`.
- Runner naming remains `home-self-host`, `home-self-host-2`, ...
- For automation/integration, emit machine-readable JSON:

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --min 1 \
  --max 4 \
  --target 4 \
  --inventory-out /tmp/runner-pool-inventory.json \
  --status-out /tmp/runner-pool-status.json
```

## 9.2) Keep 4 warm slots with always-on supervisor (recommended)

Use supervisor mode to avoid queue gaps between timer ticks.

One-command setup:

```bash
sudo -v
scripts/runner/setup_pool_supervisor_one_command.sh \
  --repo grace-bidos/chem-model-edit \
  --min 1 \
  --max 4 \
  --target 4 \
  --interval 15 \
  --token-source app \
  --app-id <APP_ID> \
  --app-installation-id <INSTALLATION_ID> \
  --app-private-key-file <PRIVATE_KEY_PATH>
```

Dry-run:

```bash
scripts/runner/setup_pool_supervisor_one_command.sh --dry-run
```

Verify:

```bash
scripts/runner/verify_self_hosted_runner_setup.sh \
  --owner grace-bidos \
  --repo chem-model-edit

sudo systemctl status chem-runner-pool-supervisor.service --no-pager
sudo systemctl status chem-github-app-token-refresh.timer --no-pager
sudo cat /var/lib/chem-model-edit/github-app-token-refresh-status.json
gh api repos/grace-bidos/chem-model-edit/actions/runners --jq '.runners[] | {name,status,busy}'
```

Optional: timer mode remains available as fallback, but supervisor mode should be the default for steady 4-slot operation.

Token lifecycle rules:

- Refresh installation tokens at least every 50 minutes for long-lived supervisors.
- Never persist app private keys in repository paths.
- Remove stale token files during incident cleanup.
- The one-command `--token-source app` path installs automated refresh by default.
- If timer setup is intentionally disabled, treat it as break-glass and set follow-up to restore automation immediately.

## 10) Emergency fallback to hosted routing

If self-hosted runners are unavailable during an incident, immediately route trusted PRs back to hosted:

- set repository variable `CI_SELF_HOSTED_TRUSTED_ROUTING=false`

This is safe to apply before or during local runner recovery.

## 11) Post-rollout success criteria and ownership

Post-rollout success criteria (evaluate on latest 50 trusted PR runs):

- runner lifecycle success rate >= 98%
- routing correctness = 100% (trusted/untrusted split)
- p95 queue wait <= 120 seconds
- p95 setup latency <= 180 seconds from job start

Ownership and escalation:

- Primary owner: repository maintainers for CI operations
- Secondary owner: infrastructure maintainer for host/VM runtime
- Escalation SLA:
  - acknowledge incident within 15 minutes
  - restore safe routing state (hosted fallback or healthy self-hosted) within 30 minutes
