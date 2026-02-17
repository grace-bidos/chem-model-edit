# Self-hosted JIT+VM Local Setup

## Goal

Run CI for trusted pull requests on self-hosted ephemeral runners backed by a local VM pool, while keeping untrusted contexts on GitHub-hosted runners.

For concrete operator commands, see:

- `docs/process/self-hosted-jit-vm-operator-checklist.md`

## Standard Ops Policy (Default)

- Primary control mode is always-on supervisor:
  - `scripts/runner/setup_pool_supervisor_one_command.sh`
- Trusted routing is enabled by default:
  - `CI_SELF_HOSTED_TRUSTED_ROUTING=true`
- Keep timer reconcile mode as fallback only:
  - `scripts/runner/setup_pool_reconcile_one_command.sh`
- Use one-command wrappers for fast operator actions:
  - token refresh automation: `scripts/runner/setup_github_app_token_refresh_one_command.sh`
  - base-runner recovery: `scripts/runner/recover_base_runner_one_command.sh`
  - routing guard rollback helper: `scripts/runner/guard_trusted_routing.sh`

## Trusted Routing Policy

Self-hosted route is selected only when all conditions match:

- event is `pull_request`
- PR is not from a fork
- actor is not a bot (`[bot]`, `dependabot[bot]`)
- `author_association` is one of `OWNER`, `MEMBER`, `COLLABORATOR`
- repo variable `CI_SELF_HOSTED_TRUSTED_ROUTING` is `true`

Fallback target is always `ubuntu-latest`.

Merge queue compatibility:

- CI must also run on `merge_group` events.
- Required checks for queue entry should remain fast (`web`, `api`, `contract`).
- Heavy suites should run post-merge on `main` to avoid slowing stacked PR throughput.

## Architecture (JIT + VM Autoscaler)

- A small controller receives queue pressure or polling signals.
- Controller boots VM from a reusable base image.
- Controller requests JIT runner config from GitHub API.
- Runner starts in ephemeral mode on the VM and executes one job.
- VM is torn down after completion.

## Warm VM Guidance

- `1 VM = 1 Job` is the safest default and is recommended.
- Warm standby (`1 idle VM`) typically saves most of boot overhead.
- For local persistent pool operation, use:
  - minimum warm concurrency: `--min` (floor)
  - max concurrency: `--max` (ceiling)
  - current desired concurrency: `--target` (optional, clamped to floor/ceiling)
  - machine-readable outputs:
    - `--inventory-out <path>` for runner inventory snapshot JSON
    - `--status-out <path>` for scaler run status JSON
  - recommended control mode: always-on supervisor
    - `scripts/runner/setup_pool_supervisor_one_command.sh`
- Practical impact target:
  - cold start overhead: around tens of seconds to low minutes depending on host and image size
  - warm standby reduces that to near runner registration time plus checkout/dependency setup

## Cache Strategy

- Keep pnpm store cache at `~/.pnpm-store`.
- Keep uv cache at `~/.cache/uv`.
- Save cache on `main` only; restore on PRs and other branches.
- Do not cache `node_modules` or `.venv`.

## Required Operator Inputs

The following requires your action:

- GitHub org/repo auth configuration for JIT issuance
- Runner group creation and repository access scope
- host sudo tasks for VM image build, service install, and lifecycle hooks

Recommended auth model for JIT issuance:

- Prefer GitHub App installation tokens (short-lived) over long-lived PATs.
- Treat token files as transient secrets (`0600`, root-owned, periodic refresh via systemd timer).
- Keep PAT only as emergency fallback during incident recovery.

## Verification Checklist

- trusted PR routes to label set `self-hosted,linux,x64,chem-trusted-pr`
- untrusted PR still runs on `ubuntu-latest`
- runner lifecycle is ephemeral (single job then removed)
- failed jobs still clean up VM and runner registration
- rollback path works by setting `CI_SELF_HOSTED_TRUSTED_ROUTING=false`
- fast rollback/recover commands are documented and tested:
  - disable: `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body false`
  - restore default: `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body true`
  - auto-guarded disable on degraded health:
    - `scripts/runner/guard_trusted_routing.sh --owner <owner> --repo <repo>`
- local health check command reports `0` when service/GitHub status align
  - `scripts/runner/check_local_runner_health.sh --owner <owner> --repo <repo>`
- local recovery command supports dry-run and reconfigure/restart flow
  - `RUNNER_OWNER=<owner> RUNNER_REPO=<repo> RUNNER_LABELS=<labels> RUNNER_GROUP=<group> scripts/runner/recover_base_runner.sh --dry-run`
  - `scripts/runner/recover_base_runner_one_command.sh --dry-run`
- local supervisor command supports dry-run before enabling
  - `scripts/runner/setup_pool_supervisor_one_command.sh --dry-run`
- app token refresh timer is active and recent status is recorded
  - `sudo systemctl status chem-github-app-token-refresh.timer --no-pager`
  - `sudo cat /var/lib/chem-model-edit/github-app-token-refresh-status.json`

## Rollback

- Immediate rollback:
  - `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body false`
- Guarded rollback (recommended operational trigger):
  - `scripts/runner/guard_trusted_routing.sh --owner <owner> --repo <repo>`
- Return to standard policy:
  - `gh variable set CI_SELF_HOSTED_TRUSTED_ROUTING --repo <owner>/<repo> --body true`
- Hard rollback: disable runner group mapping for this repository.
- Emergency fallback note: always toggle routing false first when incidents cause queueing or local runner instability.
