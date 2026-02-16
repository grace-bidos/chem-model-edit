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

## Schemathesis Playbook

- 既定の `smoke` モードは高速重視（`max-examples=8`）で、`not_a_server_error` に加えて `content_type_conformance` / `response_headers_conformance` を実行する。
- `broad` モードは探索重視（`max-examples=30`）で、`status_code_conformance` / `response_schema_conformance` / `negative_data_rejection` を含む拡張チェックセットを実行する。
- 再現性のため、既定で `--seed` と `--generation-deterministic` を有効化する（この場合 generation DB は無効）。`SCHEMATHESIS_DETERMINISTIC=0` のときは `.schemathesis/examples.db` を利用できる。
- 実行ログには mode / include regex / checks / seed / 失敗上限 / deterministic / generation DB / reports などの再実行メタデータを出力する。
- 必要に応じて `SCHEMATHESIS_REPORTS=ndjson,junit` と `SCHEMATHESIS_REPORT_DIR=...` を指定してレポートを保存する。

## Type & Quality Policy

- ローカル: できるだけ早く厳密に実行する（`api-quality-fast` + `api-security` を推奨）
- CI: 段階導入で運用する
- `pyright` / `bandit` / `pip-audit` / `gitleaks` / `deptry` / `vulture` / `mutmut` は phase-1 で non-blocking
- ノイズ削減後、段階的に required check へ昇格する
