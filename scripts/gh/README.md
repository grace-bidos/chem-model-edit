# GitHub Helper Scripts

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
