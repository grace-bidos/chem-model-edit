# GitHub Helper Scripts

## `bootstrap_pr_body.sh`

Generate a PR body stub with required metadata fields and CodeRabbit policy selection.

Examples:

```bash
scripts/gh/bootstrap_pr_body.sh GRA-60 /tmp/pr_body.md --type Ship --size XS --queue Optional --stack Standalone --coderabbit Optional
scripts/gh/create_pr_from_template.sh main feature/gra-60-pr-ci-helper-scripts "feat: add PR/CI helper scripts" /tmp/pr_body.md
```

## `pr_readiness.py`

One-shot PR readiness summary for merge checks, unresolved review threads, and branch-behind state.

Examples:

```bash
scripts/gh/pr_readiness.py 123
scripts/gh/pr_readiness.py https://github.com/<owner>/<repo>/pull/123
```

## `pr-autoloop.py`

Watch a PR, report readiness blockers, and optionally merge when ready.

Examples:

```bash
# One-shot status check
scripts/gh/pr-autoloop.py 123

# Watch and auto-merge when ready
scripts/gh/pr-autoloop.py 123 --watch --merge-when-ready --merge-method merge

# Watch with timeout and resolve only outdated threads
scripts/gh/pr-autoloop.py 123 --watch --max-wait 3600 --resolve-outdated-threads
```

## `ci_log_tail.sh`

Tail failed-step logs from a workflow run.

Examples:

```bash
scripts/gh/ci_log_tail.sh 123 --lines 200
scripts/gh/ci_log_tail.sh --run-id 123456789 --lines 80
```
