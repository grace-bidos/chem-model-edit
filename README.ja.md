# chem-model-edit

[English](./README.md) | 日本語

Quantum ESPRESSO `.in` の構造を可視化・編集するWebアプリ。

ゴール: QE `.in` 構造の view / edit / compare / share をブラウザで完結させる。  
詳細なスコープと受け入れ条件: `specs/001-chem-model-webapp/spec.md`

## クイックスタート
コマンドは repo root で実行します。

```bash
git clone git@github.com:grace-bidos/chem-model-edit.git
cd chem-model-edit

./scripts/setup-dev.sh
```

`just` は日常開発（`just dev`, `just test`, `just typecheck`）に推奨ですが、必須ではありません。

`just` の手動インストール:
```bash
# Rust/Cargo 環境が必要
cargo install just

# Linux (snap)
sudo snap install --classic just
```

`just` がない場合は個別起動:
```bash
pnpm -C apps/web dev
cd apps/api && uv run uvicorn main:app --reload --port 8000
```

セットアップ中に `cargo` 経由で `just` を自動導入したい場合:
```bash
SETUP_INSTALL_JUST=1 ./scripts/setup-dev.sh
just dev
```

- `./scripts/setup-dev.sh` は pnpm依存導入、build script承認、APIの `uv sync` を実行します。
- `just` は任意です。
- `just dev` は Web + API を同時起動します。
  - `just` 既定ポート: `WEB_PORT=3001`, `API_PORT=8000`
  - 使用中ポートがある場合は空きポートへ自動変更し、ログに表示します。

詳細セットアップ: `docs/setup.ja.md`  
ZPE worker セットアップ: `docs/zpe-worker-setup.ja.md`

## 個別起動（代替手順）
以下のコマンドは repo root で実行（特記除く）。

### Web のみ
```bash
pnpm -C apps/web dev
```
直接 Vite 起動時の既定ポートは `3000`。

### API のみ
```bash
cd apps/api
uv sync
uv run uvicorn main:app --reload --port 8000
```

### `just dev` のポート上書き
```bash
WEB_PORT=4000 API_PORT=9000 just dev
```

## 品質チェック
```bash
just style      # web: prettier + eslint, api: ruff
just typecheck  # web: tsc, api: mypy
just test       # web: vitest, api: pytest
just ci         # nx run-many -t lint,typecheck,test,knip
```

## リポジトリ構成
- `apps/web`: TanStack Start web app（frontend）
- `apps/api`: FastAPI backend
- `packages/api-client`: APIクライアント共有パッケージ
- `packages/shared`: shared types
- `docs`: セットアップ/運用ドキュメント
- `specs`: Spec/Plan/Task
- `samples/qe-in`: QE `.in` サンプル入力

## 補足（サンドボックス環境）
- Codex CLI などのサンドボックスでは、ファイルシステム制約（`EXDEV`）により `uv sync` / `just` 実行で昇格が必要な場合があります。
