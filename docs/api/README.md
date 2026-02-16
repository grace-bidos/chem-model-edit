# API Visualization

バックエンドの仕様・構造可視化は次のレシピで生成する。

- `just api-openapi-export`: OpenAPI JSON を更新
- `just api-openapi-generate`: openapi-generator で HTML/Markdown を生成
- `just api-pyreverse`: `app/` と `services/` の UML 風依存図を生成
- `just api-cov`: pytest coverage レポートを生成
- `just api-typecheck-strict`: mypy + pyright の厳密型チェック
- `just api-security`: bandit + pip-audit
- `just api-deadcode`: deptry + vulture
- `just api-schemathesis-smoke`: Schemathesis（高速スモーク）
- `just api-schemathesis-broad`: Schemathesis（広範囲探索）
- `just api-mutation-smoke`: mutmut のスモーク実行
- `just api-quality-fast`: 日常開発向けの高速品質ゲート

主な出力:

- `packages/api-client/openapi/openapi.json`
- `docs/api/openapi-html/index.html`
- `docs/api/openapi-markdown/README.md`
- `docs/graphs/api-classes.svg`
- `docs/graphs/api-packages.svg`
- `apps/api/htmlcov/index.html`
- `apps/api/coverage.xml`

## Type & Quality Policy

- ローカル: できるだけ早く厳密に実行する（`api-quality-fast` + `api-security` を推奨）
- CI: 段階導入で運用する
- `pyright` / `bandit` / `pip-audit` / `gitleaks` / `deptry` / `vulture` / `mutmut` は phase-1 で non-blocking
- ノイズ削減後、段階的に required check へ昇格する
