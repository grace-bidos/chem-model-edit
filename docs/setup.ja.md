# セットアップ（ローカル開発）

`chem-model-edit` のローカル開発環境を整える手順です。
Quantum ESPRESSO (QE) の実行環境は **不要** です。
ZPE の compute-plane（worker）セットアップは `docs/zpe-worker-setup.ja.md` を参照してください。

## 前提
- Git
- Node.js（Corepack 有効化）
- pnpm 10.27.0
- Python >= 3.13（`uv` を推奨）
- 任意: `just`（Web+API 同時起動の `just dev` を使う場合）

`just` は `just dev` / `just test` / `just typecheck` などの日常開発コマンドで推奨ですが、任意です。

`just` を入れる場合:
- Cargo 経由（Rust/Cargo 環境が必要）: `cargo install just`
- Linux の snap 経由: `sudo snap install --classic just`
- もしくは OS のパッケージマネージャ

## 推奨導線
推奨の入口は `./scripts/setup-dev.sh` です。

```bash
./scripts/setup-dev.sh
```

このセットアップで実行される内容:
- `corepack enable` + `corepack prepare pnpm@10.27.0 --activate`
- `pnpm install`
- `pnpm approve-builds`（`@parcel/watcher` を承認）
- `uv python install 3.13` と `apps/api` での `uv sync`

補足:
- `pnpm approve-builds` により `pnpm-workspace.yaml` が更新されます。
- `SETUP_SKIP_APPROVE=1` を指定すると `pnpm approve-builds` をスキップできます。
- `SETUP_INSTALL_JUST=1` を指定すると、`cargo` がある場合のみ `just` を自動導入します（Rust/Cargo 環境が必要）。

## 手動セットアップ
```bash
# Node + pnpm
corepack enable
corepack prepare pnpm@10.27.0 --activate
pnpm install

# build scripts を承認（@parcel/watcher を選択）
pnpm approve-builds

# Python + API 依存関係
uv python install 3.13
cd apps/api
uv sync
```

## 起動
### Web + API 同時起動（推奨）
```bash
just dev
```

`just dev` の実行には `just` の導入が必要です。

`just dev` のポート挙動:
- 既定値は `WEB_PORT=3001` と `API_PORT=8000`
- 既に使用中なら空きポートを自動選択してログに表示

明示的にポート指定する場合:
```bash
WEB_PORT=4000 API_PORT=9000 just dev
```

### 個別起動（代替手順）
```bash
# Web
pnpm -C apps/web dev

# API
cd apps/api
uv run uvicorn main:app --reload --port 8000
```
Web を Vite で直接起動する場合の既定ポートは `3000` です。

## 品質チェック
```bash
just style
just typecheck
just test
just ci
```

レシピの詳細は `just --list` と `Justfile` を参照してください。

## トラブルシュート
- サンドボックス環境（例: Codex CLI）では `uv sync` が `EXDEV`（cross-device rename）で失敗する場合があります。
  その場合は昇格実行で `uv sync` を行ってください。
