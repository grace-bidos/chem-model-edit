# Architecture Decision Records

This directory stores architecture decisions that are binding for backend refresh work.

- `ADR-0001`: System-of-record boundaries and ownership model (`GRA-12`)
- `ADR-0002`: Canonical JobState contract and transition rules (`GRA-13`)

Rules:

- New runtime-impacting architecture decisions must be captured in a new ADR.
- Superseding a decision requires a new ADR that references the old one.
- Implementation PRs touching these boundaries must link the relevant ADR.
