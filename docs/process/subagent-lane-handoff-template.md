# Sub-agent Lane Handoff Template

Use this template when a lane owner hands off merge-ready status to the main agent.
Keep entries concise and factual.

## Lane Summary

- Linear child issue:
- PR number and title:
- Branch / worktree:
- Lane conflict class: low / medium / high
- Scope boundary (owned files/modules):

## Merge-readiness Contract

- Required checks: pass / fail (list failing checks if any)
- Unresolved review threads: 0 / not 0
- Head branch status: not `BEHIND` / `BEHIND`
- Contract-sensitive artifacts updated in same PR: yes / no / n/a

## Review/CI Delta Since Last Handoff

- Latest review feedback addressed:
- CI failures fixed:
- Remaining blockers:

## Risk and Conflict Notes

- Cross-lane conflict detected: yes / no
- If yes, impacted files:
- Proposed resolution options:

## Linear and Queue State

- Linear state: `In Review` / `Done` (with acceptance criteria note)
- Queue policy in PR body: Required / Optional
- Merge queue action taken:

## Post-merge Cleanup Plan

- Worktree cleanup command:
- Branch cleanup command:
- Local main sync command:
