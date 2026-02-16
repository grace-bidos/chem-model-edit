# VM Validation Lane Secret Rotation Playbook (GRA-123)

This playbook defines what to do when VM validation lane credentials are exposed, suspected exposed, or intentionally rotated.

## Scope

- VM validation lanes and helper scripts used by lane owners.
- Local temporary artifacts for JIT config and VM credentials.
- Process-only guidance; do not store real secrets in the repository.

## Trigger Conditions

Run this playbook when any of the following occurs:

- secret material appears in terminal history, logs, screenshots, or PR comments
- temporary credential files were created in the repository by mistake
- a key/token age policy requires scheduled rotation
- credentials are shared across operators and ownership changes

## Immediate Response (First 15 Minutes)

1. Stop active automation that may still consume leaked credentials.
2. Revoke or rotate exposed credentials in the source system first.
3. Invalidate all local temporary files and shell variables.
4. Replace local files with newly issued credentials only in local ignored paths.
5. Record the incident in the linked Linear child issue and describe blast radius.

## Safe Local Handling Rules

- Store temporary VM/JIT credentials only under:
  - `.tmp/vm-validation-secrets/`
  - `.tmp/vm-jit-config/`
  - `.secrets/vm-validation/`
- Never commit temporary credentials, even if encrypted.
- Do not paste raw secret values into PR descriptions, issue comments, or CI logs.
- Prefer short-lived credentials and one-time artifacts.
- Remove credentials after each validation run.

## Rotation Checklist

Use this checklist for every rotation event:

- [ ] Old credential revoked/disabled in provider control plane.
- [ ] New credential created with least privilege.
- [ ] Runtime consumer updated (VM, runner, or operator environment).
- [ ] Validation lane rerun confirms healthy authentication.
- [ ] Local temporary files cleaned with `scripts/runner/cleanup_vm_temp_credentials.sh`.
- [ ] Incident notes added to the related Linear issue (no raw secret values).

## Verification Commands

Basic local verification before/after rotation:

```bash
# Dry-run cleanup visibility
scripts/runner/cleanup_vm_temp_credentials.sh --dry-run

# Perform cleanup without prompt (CI-safe/local automation)
scripts/runner/cleanup_vm_temp_credentials.sh --force
```

## Local Artifact Cleanup Policy

- Use the cleanup script after each VM lane run and after any failed attempt.
- If manual cleanup is required, delete only files under ignored local paths.
- If an artifact location is outside allowed local paths, stop and handle manually.

## Reporting Requirements

When rotation happens, include in the lane handoff:

- what was rotated (token/key category only)
- when it was rotated (UTC timestamp)
- what systems were updated
- confirmation that local temporary artifacts were removed

Do not include secret values, partial values, or screenshots containing credentials.
