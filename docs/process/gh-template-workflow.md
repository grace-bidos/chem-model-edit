# GitHub Template Workflow

Use body files instead of inline shell strings for issue/PR creation to avoid shell expansion bugs.

## Create Issue

```bash
scripts/gh/create_issue_from_template.sh "Title" /tmp/issue_body.md task
```

Use this only when creating an optional GitHub mirror/discussion issue.
Linear remains planning source-of-truth.

## Create PR

```bash
scripts/gh/create_pr_from_template.sh <base-branch-or-main> <head-branch> "PR title" /tmp/pr_body.md --draft
```

## Notes

- Keep public issue/PR text in English.
- Prefer markdown files committed in `specs/` or `/tmp/*.md` when drafting bodies.
- Ensure PR body includes required metadata fields checked by PR policy (Linear issue, type, size, queue policy, stack, and CodeRabbit policy).
