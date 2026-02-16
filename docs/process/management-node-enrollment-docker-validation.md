# Management Node Enrollment Docker Validation (GRA-82)

This document records the lane-B validation flow for management-node enrollment in Docker, using bootstrap + Ansible artifacts recreated under this lane.

## Scope

- Bootstrap entrypoint: `scripts/bootstrap-management-node.sh`
- Ansible project root: `ops/ansible/management-node-enrollment`
- Docker harness root: `ops/ansible/management-node-enrollment/docker-harness`
- Validation script: `scripts/verify-management-node-enrollment-docker.sh`

## Docker Harness

Run from repository root:

```bash
./scripts/verify-management-node-enrollment-docker.sh
```

The script performs:

1. Harness reset (`docker compose down --volumes --remove-orphans`)
2. Harness build/start (`docker compose up -d --build`)
3. In-container tool checks (`ansible --version`, `python --version`, `bash -n`)
4. Enrollment dry-run (`bootstrap-management-node.sh` default mode)
5. Enrollment apply (`bootstrap-management-node.sh --apply`)
6. Enrollment idempotency apply (`bootstrap-management-node.sh --apply` again)

## Success Criteria

- Dry-run exits `0`
- Apply exits `0`
- Idempotency apply exits `0`
- Idempotency recap contains `changed=0 unreachable=0 failed=0`

## Notes

- The role intentionally skips `authorized_key` in check mode by default via `management_enrollment_skip_key_install_in_check_mode: true` to avoid a known check-mode ordering limitation when user creation is also simulated.
- Apply mode still performs key installation and sudoers updates.
- Harness image includes `sudo` so sudoers validation (`visudo -cf`) can run in-container.
