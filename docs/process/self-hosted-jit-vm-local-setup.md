# Self-hosted JIT+VM Local Setup

## Goal

Run CI for trusted pull requests on self-hosted ephemeral runners backed by a local VM pool, while keeping untrusted contexts on GitHub-hosted runners.

## Rollout Model

- Phase 0: keep `CI_SELF_HOSTED_TRUSTED_ROUTING=false` (no traffic, routing logic present).
- Phase 1: set `CI_SELF_HOSTED_TRUSTED_ROUTING=true` and canary only trusted PRs.
- Phase 2: tune warm VM policy and cache behavior.

## Trusted Routing Policy

Self-hosted route is selected only when all conditions match:

- event is `pull_request`
- PR is not from a fork
- actor is not a bot (`[bot]`, `dependabot[bot]`)
- `author_association` is one of `OWNER`, `MEMBER`, `COLLABORATOR`
- repo variable `CI_SELF_HOSTED_TRUSTED_ROUTING` is `true`

Fallback target is always `ubuntu-latest`.

## Architecture (JIT + VM Autoscaler)

- A small controller receives queue pressure or polling signals.
- Controller boots VM from a reusable base image.
- Controller requests JIT runner config from GitHub API.
- Runner starts in ephemeral mode on the VM and executes one job.
- VM is torn down after completion.

## Warm VM Guidance

- `1 VM = 1 Job` is the safest default and is recommended.
- Warm standby (`1 idle VM`) typically saves most of boot overhead.
- Practical impact target:
  - cold start overhead: around tens of seconds to low minutes depending on host and image size
  - warm standby reduces that to near runner registration time plus checkout/dependency setup

## Cache Strategy

- Keep pnpm cache on `actions/setup-node`.
- Enable uv cache on `astral-sh/setup-uv`.
- Save cache primarily on `main`; restore on PRs.
- Do not cache `node_modules` or `.venv`.

## Required Operator Inputs

The following requires your action:

- GitHub org/repo auth configuration for JIT issuance
- Runner group creation and repository access scope
- host sudo tasks for VM image build, service install, and lifecycle hooks

## Verification Checklist

- trusted PR routes to label set `self-hosted,linux,x64,chem-trusted-pr`
- untrusted PR still runs on `ubuntu-latest`
- runner lifecycle is ephemeral (single job then removed)
- failed jobs still clean up VM and runner registration
- rollback path works by setting `CI_SELF_HOSTED_TRUSTED_ROUTING=false`

## Rollback

- Immediate rollback: set `CI_SELF_HOSTED_TRUSTED_ROUTING=false`.
- Hard rollback: disable runner group mapping for this repository.
