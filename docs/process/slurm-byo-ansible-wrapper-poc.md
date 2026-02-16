# Slurm BYO Ansible Wrapper PoC (GRA-71)

This PoC adds a safe-by-default wrapper that converts onboarding inputs into Ansible inventory and vars files for a minimal control+worker Slurm topology.

## Scope

- Wrapper script: `scripts/bootstrap-slurm-byo.sh`
- Sample Ansible harness: `ops/ansible/slurm-byo/examples`
- Output artifacts:
  - inventory with `slurmservers` and `slurmexechosts`
  - vars shaped for `galaxyproject.slurm` baseline (`slurm_roles`, `slurm_config`, `slurm_nodes`, `slurm_partitions`)

## Baseline Reference (galaxyproject/ansible-slurm)

The PoC follows baseline role conventions from `galaxyproject.slurm`:

- Inventory groups: `slurmservers`, `slurmexechosts`.
- Role-selection variable: `slurm_roles` (`controller` for control plane, `exec` for workers when expanded).
- Core config structures: `slurm_config`, `slurm_nodes`, and `slurm_partitions`.

This PoC intentionally keeps execution tasks to `debug` output in sample playbooks. Real cluster bootstrap should replace those tasks with the actual role invocation and pinned role versioning.

## Usage

Plan-only (generate files, print command, do not run ansible):

```bash
scripts/bootstrap-slurm-byo.sh \
  --plan-only \
  --cluster-name chem-byo \
  --controller-host 10.0.0.10 \
  --worker-hosts 10.0.0.11,10.0.0.12 \
  --ssh-user ubuntu
```

Dry-run check mode (default if `--apply` is not set):

```bash
scripts/bootstrap-slurm-byo.sh \
  --cluster-name chem-byo \
  --controller-host 10.0.0.10 \
  --worker-host 10.0.0.11 \
  --worker-host 10.0.0.12 \
  --ssh-user ubuntu
```

Apply mode (opt-in only):

```bash
scripts/bootstrap-slurm-byo.sh \
  --apply \
  --cluster-name chem-byo \
  --controller-host 10.0.0.10 \
  --worker-hosts 10.0.0.11,10.0.0.12 \
  --ssh-user ubuntu
```

## Idempotent Rerun Flow

1. Run with `--plan-only` to verify generated inventory/vars and command line.
2. Run default dry-run mode and inspect check-mode output.
3. If acceptable, run with `--apply`.
4. Rerun with the same inputs after apply. The generated files are deterministic and replaceable, so re-execution should converge.
5. For topology changes (add/remove workers), rerun with updated inputs and inspect diff in generated vars before apply.

## Risks and Limits (Advanced/HA)

- HA controllers are not modeled: no multi-controller quorum, failover VIP, or split-brain protections.
- No `slurmdbd`/accounting DB orchestration is included.
- No OS/package prerequisites are enforced in wrapper-generated artifacts.
- Worker aliases are generated as `workerN`; advanced host naming or mixed pools may need template extension.
- Secrets handling is intentionally minimal for PoC; production use should move credentials/keys to Vault or equivalent.
- Inventory and vars generation is local-file based; concurrent operators can overwrite shared runtime paths unless path overrides are used.

## Promotion Path

- Add pinned `galaxyproject.slurm` role dependency and role execution playbook.
- Split controller/worker role vars cleanly (including `slurm_roles: ['exec']` on workers).
- Add HA-aware inventory model and explicit validation for controller count/topology.
- Introduce CI checks that run wrapper in plan-only and check-mode paths with fixture inputs.
