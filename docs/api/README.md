# API Visualization

## Testing

API tests are executed with `pytest` profiles below.

- `just api-test`: default full API test run (backward-compatible command)
- `just api-test-fast`: local fast profile (`-n auto`, excludes `e2e/contract/slow/schemathesis`)
- `just api-test-coverage`: coverage profile, generates HTML/XML reports (`--cov-fail-under=75.5`, phase-1 gate)
- `just api-cov`: compatibility alias to `just api-test-coverage`
- `pnpm nx test api --configuration=fast`: Nx fast profile
- `pnpm nx test api --configuration=coverage`: Nx coverage profile

Coverage gate roadmap (phase plan):

- Phase 1 (current): `--cov-fail-under=75.5` (locked to current measured baseline)
- Phase 2 (next cycle target): `--cov-fail-under=76`
- Phase 3 (target): `--cov-fail-under=80`

Plugin note:

- `pytest-time-machine` package does not exist on PyPI; this repo uses `time-machine` (pytest-compatible) for time-travel tests.

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

## Schemathesis (`/api/*` 全パス対象)

`apps/api/scripts/schemathesis.sh` は OpenAPI の `/api/*` 全パスを対象に実行する。
デフォルトは CI 向けの軽量 `smoke` モードで、`SCHEMATHESIS_MODE=deep` を指定すると探索量を増やせる。

基本実行:

- `pnpm exec nx run api:schemathesis`

主な環境変数:

- `SCHEMATHESIS_MODE=smoke|deep` (`smoke` 既定)
- `SCHEMATHESIS_SEED=<int>` (既定: `137`, 再現実行用)
- `SCHEMATHESIS_MAX_EXAMPLES=<int>` / `SCHEMATHESIS_MAX_FAILURES=<int>` (モード既定の上書き)
- `SCHEMATHESIS_AUTH_TOKEN=<token>` または `SCHEMATHESIS_AUTH_HEADER='Authorization: Bearer <token>'`
- `SCHEMATHESIS_TENANT_ID=<tenant-id>` (`x-tenant-id` ヘッダーとして付与)
- `SCHEMATHESIS_DEV_USER_ID=<id>` / `SCHEMATHESIS_DEV_USER_EMAIL=<email>` (dev-bypass ヘッダー)

認証付きエンドポイントをローカル検証する例:

- `AUTH_MODE=dev-bypass ZPE_ADMIN_TOKEN=secret SCHEMATHESIS_DEV_USER_ID=user-1 SCHEMATHESIS_DEV_USER_EMAIL=user-1@example.com SCHEMATHESIS_TENANT_ID=tenant-dev SCHEMATHESIS_AUTH_TOKEN=secret pnpm exec nx run api:schemathesis`

深掘り実行の例:

- `SCHEMATHESIS_MODE=deep SCHEMATHESIS_SEED=137 SCHEMATHESIS_MAX_EXAMPLES=50 pnpm exec nx run api:schemathesis`
