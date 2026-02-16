# GRA-71 Investigation: Ansible wrapper PoC for BYO Slurm bootstrap

## Objective

Validate a wrapper flow that transforms onboarding inputs into Ansible inventory and vars for a minimal Slurm control+worker setup, with safe defaults.

## Implemented Artifacts

- `scripts/bootstrap-slurm-byo.sh`
- `ops/ansible/slurm-byo/examples/ansible.cfg`
- `ops/ansible/slurm-byo/examples/playbooks/site.yml`
- `ops/ansible/slurm-byo/examples/inventory/hosts.example.ini`
- `ops/ansible/slurm-byo/examples/group_vars/all.example.yml`
- `docs/process/slurm-byo-ansible-wrapper-poc.md`

## Validation Commands and Summary

### 1) Shell syntax check

Command:

```bash
bash -n scripts/bootstrap-slurm-byo.sh
```

Result:

- Pass (no syntax errors).

### 2) Plan-only path

Command:

```bash
scripts/bootstrap-slurm-byo.sh \
  --plan-only \
  --cluster-name chem-byo-poc \
  --controller-host 127.0.0.1 \
  --worker-hosts 127.0.0.1,127.0.0.1 \
  --ssh-user "$USER"
```

Observed summary:

- Wrapper generated inventory and vars into `.just-runtime/slurm-byo-ansible-wrapper/`.
- Wrapper printed a complete `ansible-playbook` command.
- No ansible execution occurred (as expected for plan-only).

### 3) Dry-run/check-mode path

Command:

```bash
scripts/bootstrap-slurm-byo.sh \
  --cluster-name chem-byo-poc \
  --controller-host 127.0.0.1 \
  --worker-hosts 127.0.0.1,127.0.0.1 \
  --ssh-user "$USER"
```

Observed summary:

- Wrapper executed ansible with `--check --diff`.
- Sample playbook completed and printed mapped Slurm baseline vars.
- No apply changes were requested because dry-run is default.

## Notes

- This PoC is intentionally limited to a minimal control+worker topology and generated artifacts.
- HA, controller failover, and accounting DB deployment are documented as explicit non-goals/risk areas.
