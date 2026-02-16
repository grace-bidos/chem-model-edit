# GRA-82 Investigation: Management-node Enrollment Validation in Docker

## Summary

Validated management-node enrollment flow end-to-end in Docker with:

- dry-run
- apply
- idempotency apply

Final verification run completed successfully on **2026-02-16**.

## Artifacts Recreated in This Lane

- `scripts/bootstrap-management-node.sh`
- `scripts/verify-management-node-enrollment-docker.sh`
- `ops/ansible/management-node-enrollment/*`
- `ops/ansible/management-node-enrollment/docker-harness/*`

## Failures Encountered and Fixes

### 1) Role resolution failure

Error:

- `the role 'management_node_enrollment' was not found`

Root cause:

- `ansible-playbook` did not consistently use the project-local `ansible.cfg` (missing `roles_path` context).

Fix:

- Updated bootstrap command to force `ANSIBLE_CONFIG`:
  - `env ANSIBLE_CONFIG=$ANSIBLE_DIR/ansible.cfg ansible-playbook ...`

### 2) Callback plugin incompatibility

Error:

- `The 'community.general.yaml' callback plugin has been removed`

Root cause:

- Legacy config used `stdout_callback = yaml`, which is not supported in the current Ansible core.

Fix:

- Updated configs to:
  - `stdout_callback = default`
  - `result_format = yaml`

### 3) localhost SSH unreachable in harness

Error:

- `Failed to connect to the host via ssh: ... localhost port 22: Connection refused`

Root cause:

- Generated inventory defaulted to SSH transport even for localhost.

Fix:

- In bootstrap-generated inventory, when host is `localhost` / `127.0.0.1` / `::1`, set:
  - `ansible_connection=local`

### 4) check-mode failure on authorized_key

Error:

- `Either user must exist or you must provide full path to key file in check mode`

Root cause:

- In dry-run, user creation is simulated and `authorized_key` cannot proceed against a not-yet-created user.

Fix:

- Added default variable:
  - `management_enrollment_skip_key_install_in_check_mode: true`
- Guarded key task with:
  - `not (ansible_check_mode and (management_enrollment_skip_key_install_in_check_mode | bool))`

## Verification Logs

From successful run:

- dry-run: `.just-runtime/management-node-enrollment/verification-logs/20260216-090852-dry-run.log`
- apply: `.just-runtime/management-node-enrollment/verification-logs/20260216-090852-apply.log`
- idempotency: `.just-runtime/management-node-enrollment/verification-logs/20260216-090852-idempotency.log`

Key recap from idempotency run:

- `changed=0`
- `unreachable=0`
- `failed=0`

## Remaining Notes

- Ansible emitted a deprecation warning from `authorized_key` internals (`to_native` import path). This is upstream module behavior and did not affect apply/idempotency correctness.
