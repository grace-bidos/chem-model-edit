# Feature Specification: 開発基盤整備 (Nx/Storybook/CI/CD/Docstring/Refactor)

> [!NOTE]
> Historical spec document. Superseded by the Modal-first runtime/deploy operating model on February 18, 2026.

**Feature Branch**: `[stack/nx-roadmap]`  
**Created**: 2026-01-31  
**Status**: Draft  
**Input**: ユーザー要望「Nx導入→Storybook→Just編集→CI(Axe-core/Fast-check)→CD→docstring/TSDoc→リファクタリング」

## 背景 / 目的

- 既存の開発フローに対して，Nxなどの導入とテストツールなどの導入を行ったあとに、統一コマンドとCI/CDの整備を行い，堅牢性を高めること。
- TanStack Startの設定が十分活用されていないため、段階的な見直しとUI基盤の整理を進めたい。

## スコープ

- Nx導入（pnpm前提、apps/packagesのproject化、相対パス運用）
- Vitest + React Testing Library の導入とNxタスク統合
- Knip導入とNxタスク統合
- Storybook導入（Chromatic連携を含む）
- Justfile更新（新ツール群向けの統一コマンド追加）
- CI整備（GitHub Actions、Axe-core/Fast-check/Schemathesis含む）
- CDの既存連携の維持（historical: Web: Cloudflare Pages、API: Cloud Run）
- docstring/TSDoc方針の整備と主要領域への適用
- リファクタリング（TanStack Start設定再構築→Shadcn table統合→その他）

## 非ゴール

- Cloudflare Pages / Cloud Run のホスティング方式の変更
- UIデザインの全面刷新（必要な範囲の改善のみ）
- API仕様の抜本変更
- 既存機能の挙動変更を伴うリライト

## 方針 / 決定事項

- CIはGitHub Actions。
- CDは既存連携を維持（Web: Cloudflare Pages、API: Cloud Run）。
- StorybookはChromaticで公開・レビューを行う。
- Axe-coreはWebのUI/コンポーネントに限定。
- Fast-checkはWeb + sharedのJS/TSロジックに限定。
- APIのプロパティテストはPython側でSchemathesisを導入。
- pnpmのパスは相対（`pnpm -C`/ワークスペース相対）で統一。

## 成果物（期待される状態）

1. **Nx基盤**
   - `nx.json` と各projectの `project.json` が配置される。
   - `apps/web`, `apps/api`, `packages/shared` がNxプロジェクトとして定義される。
   - Nxタスクは既存コマンドを呼び出す形で統一（pnpm/uvは相対パス）。
2. **Testing基盤**
   - `apps/web` にVitest + RTLが導入される。
   - A11yテストはAxe-coreベースで実装（Vitest内で実行）。
   - Fast-checkのプロパティテストを `packages/shared` もしくは `apps/web` に配置。
   - `apps/api` にSchemathesisが導入され、OpenAPIに基づいたテストが実行できる。
3. **Storybook + Chromatic**
   - `apps/web` にStorybook構成が追加される。
   - Chromatic実行用のタスク/CIステップが追加される（トークンはSecrets前提）。
4. **Justfile**
   - Nx/Storybook/Testing/Knip/CI向けのタスクが追加される。
   - 既存コマンドとの整合性が保たれ、相対パスで実行できる。
5. **CI**
   - GitHub Actionsでlint/typecheck/test/knip/storybook/axe/fast-check/schemathesisが実行される。
   - Nxを用いたキャッシュや影響範囲実行が設定される（可能な範囲）。
6. **CD**
   - 既存のCloudflare Pages / Cloud Runの連携が壊れていないことを確認する。
   - CI側の変更がCDに影響しない。
7. **Docstring/TSDoc**
   - TS/JS側にTSDoc方針を整理したドキュメントを追加。
   - 公開APIまたは共有ユーティリティにTSDocを付与。
8. **リファクタリング**
   - TanStack Start設定が現在の利用実態に合う形で整理される。
   - Shadcn tableがTable系UIで適切に統合される。
   - その他の整理は必要最小限の差分で段階的に実施。

## 成功条件

- Nxコマンドで主要タスクが一貫して実行できる。
- CIでAxe-core / Fast-check / Schemathesisが動作し、失敗時に原因が特定できる。
- StorybookがChromaticで公開・レビューできる。
- 既存のCD連携が維持される。

## 受け入れ基準

- `nx graph` で `apps/web`, `apps/api`, `packages/shared` が認識される。
- `nx run web:test`（Vitest + RTL）が成功する。
- `nx run web:a11y` でAxe-coreテストが実行される。
- `nx run shared:test` でFast-checkテストが実行される。
- `nx run api:schemathesis` がローカル/CIで実行可能。
- `nx run web:storybook` と `nx run web:chromatic` が実行可能。
- `just` タスクで上記を短縮実行できる。
- CIがGitHub Actions上で通る（Secrets/環境変数を除く）。

## 検証方法

- `nx run-many -t lint,typecheck,test,knip` を実行し通過すること。
- `nx run web:storybook` でStorybook起動を確認。
- `nx run web:chromatic --dry-run` で実行確認（本番はSecrets次第）。
- `nx run api:schemathesis` がローカルで実行できること。
- CIワークフローで同等のコマンドが通ること。

## 依存関係 / 前提

- GitHub ActionsのSecretsに `CHROMATIC_PROJECT_TOKEN` を設定する。
- Cloudflare Pages / Cloud Runの既存連携は現状維持とする。
- OpenAPIがFastAPIで有効化されていること（Schemathesis用）。

## 未確定事項

- Nxのプラグイン構成（@nx/js/@nx/vite/@nx/react など）は最小構成で合意が必要。
- Schemathesisの対象API範囲（全エンドポイント or 主要APIのみ）。
- Axe-coreの対象コンポーネント範囲（最低限の代表UI）。
