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

## 5) Tear-down guarantees

- Enforce `1 VM = 1 Job`.
- Destroy VM immediately after runner exits.
- Do not reuse runner workspace between jobs.

## 6) Enable traffic (canary)

After one successful dry run:

- set `CI_SELF_HOSTED_TRUSTED_ROUTING=true`
- open a trusted PR and verify route job selects self-hosted label target

## 6.1) Scale local runner count quickly

Use this helper to converge local runner instances to a target count:

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --target 2
```

Increase parallelism:

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --target 4
```

Dry-run:

```bash
scripts/runner/scale_local_runner_pool.sh \
  --repo grace-bidos/chem-model-edit \
  --target 4 \
  --dry-run
```

## 7) Fast rollback

Immediate rollback:

- set `CI_SELF_HOSTED_TRUSTED_ROUTING=false`

Hard rollback:

- remove repository access from the dedicated runner group.
