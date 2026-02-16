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
