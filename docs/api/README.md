# API Visualization

## Testing

API tests are executed with `pytest` profiles below.

- `just api-test`: default full API test run (backward-compatible command)
- `just api-test-fast`: local fast profile (`-n auto`, excludes `e2e/contract/slow/schemathesis`)
- `just api-test-coverage`: coverage profile, generates HTML/XML reports
- `just api-cov`: compatibility alias to `just api-test-coverage`
- `pnpm nx test api --configuration=fast`: Nx fast profile
- `pnpm nx test api --configuration=coverage`: Nx coverage profile

バックエンドの仕様・構造可視化は次のレシピで生成する。

- `just api-openapi-export`: OpenAPI JSON を更新
- `just api-openapi-generate`: openapi-generator で HTML/Markdown を生成
- `just api-pyreverse`: `app/` と `services/` の UML 風依存図を生成
- `just api-cov`: pytest coverage レポートを生成

主な出力:

- `packages/api-client/openapi/openapi.json`
- `docs/api/openapi-html/index.html`
- `docs/api/openapi-markdown/README.md`
- `docs/graphs/api-classes.svg`
- `docs/graphs/api-packages.svg`
- `apps/api/htmlcov/index.html`
- `apps/api/coverage.xml`
