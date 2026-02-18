# Local Dev Runtime Ops (WSL + Hyper-V + Modal)

This runbook is a local-operator guide for stable day-to-day development.
It keeps VM access and runtime env handling in one place, without changing runtime contracts.

## Goal

- Keep local dev secrets in one canonical file.
- Avoid repeated manual edits across `apps/api/.env` and ad-hoc files.
- Reduce WSL/Hyper-V networking confusion before VM-related tasks.

## Canonical Local Files

- Secrets source of truth (local-only): `.tmp/dev-secrets.env`
- API runtime env target: `apps/api/.env`
- Optional operator cache files: `.tmp/*`

Do not commit `.tmp/dev-secrets.env`.

## One-Time Bootstrap

```bash
scripts/dev/bootstrap_dev_secrets_file.sh
```

Then fill values in `.tmp/dev-secrets.env`.

Minimum keys for runtime submit smoke:

- `MODAL_API_BASE_URL`
- `MODAL_PROXY_KEY`
- `MODAL_PROXY_SECRET`
- `CLERK_TEST_JWT`

Additional keys for runtime gateway bridging:

- `RUNTIME_COMMAND_SUBMIT_URL`
- `RUNTIME_SERVICE_AUTH_BEARER_TOKEN`

VM lane convenience keys:

- `REMOTE_VM_HOST` (example: `grace@172.17.112.47`)
- `REMOTE_VM_REPO_DIR` (example: `/home/grace/projects/chem-model-edit`)
- `REMOTE_VM_SSH_KEY` (example: `/home/grace/.ssh/chem_baseline_ed25519`)

## Apply Runtime Settings to API Env

```bash
scripts/dev/apply_dev_runtime_env.sh
```

This updates only:

- `RUNTIME_COMMAND_SUBMIT_URL`
- `RUNTIME_SERVICE_AUTH_BEARER_TOKEN`

## Submit Smoke (FastAPI on Modal)

```bash
scripts/dev/runtime_submit_smoke.sh
```

Notes:

- `requested_by.user_id` is auto-derived from JWT `sub` in this script.
- If JWT expired, refresh `CLERK_TEST_JWT` and rerun.

## WSL + Hyper-V Troubleshooting

### Symptom: `No route to host` to VM IP in `172.17.x.x`

Cause: WSL Docker bridge often owns `172.17.0.0/16` and steals routing.

Quick check:

```bash
ip route get <vm-ip>
```

If it resolves via `docker0`, use one of these paths:

- Preferred: use Windows OpenSSH (`C:\Windows\System32\OpenSSH\ssh.exe`) in VM scripts.
- Alternative: add host route in WSL (requires sudo), then retry ssh.

### Find VM and candidate IP on Windows host (from WSL)

```bash
powershell.exe -NoProfile -Command "Get-VM | Select-Object Name,State | Format-Table -AutoSize"
powershell.exe -NoProfile -Command "arp -a"
```

For Hyper-V `Default Switch`, ARP is often the easiest way to map VM MAC to IP.

## Recommended Daily Flow

1. `scripts/dev/bootstrap_dev_secrets_file.sh` (first time or when keys changed)
2. Update `.tmp/dev-secrets.env`
3. `scripts/dev/apply_dev_runtime_env.sh`
4. Deploy API if required (`uv run --project apps/api modal deploy apps/api/modal_app.py -e dev`)
5. `scripts/dev/runtime_submit_smoke.sh`

## Security Hygiene

- Never paste raw secrets in issue/PR comments.
- Keep secrets in `.tmp/dev-secrets.env` only.
- Rotate `RUNTIME_SERVICE_AUTH_BEARER_TOKEN` when sharing environments.
