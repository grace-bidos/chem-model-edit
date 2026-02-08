# Codex Agent Context (chem-model-edit)

## 目的 / 概要
計算化学向けの構造可視化・編集ツールをWebアプリとして提供する。複数構造の並列表示、部分的な構造移植、比較・整列、エクスポート共有を支援する。

## スコープ（現時点）
- Quantum ESPRESSOの`.in`を対象に構造を読み取り・編集・書き出し
- 初期は座標と原子種のみを扱う
- 共有は「単一HTML」(NGL埋め込み) を提供。将来的にサーバ保存リンクを優先
- 3D表示は Mol* を使い、**デフォルトはボール＆スティック**

## 技術スタック
- Frontend: TanStack Start (SPA) + shadcn/ui + Mol*
- Backend: FastAPI + ASE / pymatgen
- JS/TS: pnpm
- Python: uv (+ ruff, mypy, pytest)
- Tooling: Nx, Storybook, Chromatic

## デプロイ
- Web: Cloudflare Workers
- API: Cloud Run

## リポジトリ構成
- `apps/web`: Webアプリ (TanStack Start)
- `apps/api`: FastAPI API
- `packages/shared`: 共通型
- `specs/001-chem-model-webapp`: Spec/Plan/Task
- `ref-legacy`: 旧実装 (参照用)
- `samples/qe-in`: テスト用の Quantum ESPRESSO `.in` サンプル
- `uploads`: ユーザーが外部ドキュメント/データ/サンプルをアップロードする場所

## 重要なUI/機能 (web)
- 複数構造管理、選択コピー/貼り付け、距離表示、シフト/整列
- Import: ファイル/クリップボード
- Export: QE `.in` / Share HTML
- Supercell生成: API呼び出し + Mol*プレビュー

## エージェント作業方針
- **日本語で報告**
- GitHub の PR / Issue / コメント等、外部公開される成果物は英語で記載する
- あらゆる作業の開始前に Git worktree を作成し、環境を分離してから作業する（並列作業のため）
- Spec Kit ワークフロー順守:
  - Spec → Plan → Task を先に整備
  - 変更が想定外ならSpecを更新して再承認
- 既存の仕様/タスクに従う
- `uploads` 配下のユーザーアップロード済みデータは、作業上必要であれば移動・整理・削除してよい

## 開発/検証コマンド
- Web:
  - `pnpm -C apps/web dev --port 3001`
  - `just web`（Codexサンドボックスでは昇格実行が必要）
  - `pnpm -C apps/web typecheck`
  - `pnpm -C apps/web storybook`
  - `pnpm -C apps/web chromatic`
- API:
  - `just deps`（初回/依存更新時）
  - `uv run uvicorn apps.api.main:app --reload --port 8000`
  - `uv run pytest`
  - `uv run mypy .`
- Nx:
  - `pnpm exec nx graph`
  - `pnpm exec nx run-many -t lint,typecheck,test,knip`
- Just:
  - `just style`（web: prettier+eslint / api: ruff）
  - `just test`（web: vitest / api: pytest）
  - `just typecheck`（web: tsc / api: mypy）
  - `just storybook`
  - `just chromatic`

## 既知の注意点
- Mol*は`Viewer.create`を使い、PDB読み込み後に`ball-and-stick`表現を追加
- Import/Exportの入出力は`.in` (QE) と座標/原子種に限定
- `apps/web/src/components/molstar/MolstarViewer.tsx` が表示の中心
- pnpmストアは `.pnpm-store` を利用する方針（相対パス／`.pnpm-store` はgit ignore）
- Codexサンドボックス（workspace-write）では `uv sync` / `just` 実行時は昇格（on-request）が必要

## ブランチ運用
- トランクは `main` とし、作業ブランチは必ず `main` から **短命で** 切る
- 作業ブランチごとに `git worktree` を作成し、並行作業を分離する
- `dev` のような長命ブランチは持たない（必要なら一時的な検証用途に限定）
- 基本は短期ブランチ → PR → **merge commit** で `main` に統合（squash/rebase は使わない）
- `main` への直接pushは禁止（PR経由で統合）
- レビューは必須ではない（自己マージ可）
- CIチェックは必須ではない（任意）
- 命名規約: `feat/`, `fix/`, `refactor/`, `chore/`, `docs/`, `spike/`, `codex/`, `ui/`, `stack/`
- マージ後はブランチを削除（GitHubの自動削除を有効）

## worktree運用
- worktreeの置き場は `chem-model-edit/.worktrees/<name>` に固定する（並列ディレクトリを作らない）
- `main` の作業ツリーはレビュー/確認専用に保ち、実作業は worktree で行う
- 原則「ブランチ名 = worktree名」。同一ブランチで複数必要なら `<branch>@<purpose>` を使う
- PRがマージされたら対応worktreeは削除する（`git worktree remove` → `git worktree prune`）
- ルートの `.gitignore` に `.worktrees/` を含める

## 進め方のベストプラクティス
- 変更前に影響範囲と対象ファイルを明確化
- 破壊的変更は避け、既存実装との互換性を維持
- 必要に応じてテスト/型チェックを実行
- 実装結果はSpec/Plan/Taskへ反映
- Chrome-devtools-mcpを使ってブラウザでの動作検証
- 作業の区切りのたびにgitを利用してください
- 編集後は型チェック，Lint, テストを行って

## APIレイヤ命名規約
- 層はディレクトリで区別する（例: `app/routers`, `services`）
- 同一ドメイン名は層をまたいで同じファイル名を使う（例: `routers/supercells.py` と `services/supercells.py`）
- 単数形/複数形の揺れで区別しない。判別はディレクトリ責務で行う
