# AGENTS

このリポジトリのエージェント運用方針は `Agent.md` を参照。

## APIレイヤ命名規約
- 層はディレクトリで区別する（例: `app/routers`, `services`）
- 同一ドメイン名は層をまたいで同じファイル名を使う（例: `routers/supercells.py` と `services/supercells.py`）
- 単数形/複数形の揺れで区別しない。判別はディレクトリ責務で行う
