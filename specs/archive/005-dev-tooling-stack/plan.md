# Plan: 開発基盤整備 (Nx/Storybook/CI/CD/Docstring/Refactor)

> [!NOTE]
> Historical plan document. Superseded by the Modal-first runtime/deploy operating model on February 18, 2026.

## 目的

開発基盤を統一し、CI/CD・品質・開発体験を段階的に強化する。

## 作業ステップ

1. Specレビュー（ヒューマン承認）。
2. Nx導入（pnpm前提、project化、基本ターゲット整備）。
3. Vitest/RTL + Knip導入とNx統合。
4. Storybook導入 + Chromatic連携。
5. Justfile整理（新ツール群の統一コマンド）。
6. CI整備（Axe-core/Fast-check/Schemathesis含む）。
7. CD既存連携の確認（historical: Cloudflare Pages / Cloud Run）。
8. docstring/TSDoc方針整備と適用。
9. リファクタリング（TanStack Start設定→Shadcn table→その他）。

## 完了条件

- Specが承認されている。
- 受け入れ基準を満たす。

## 依存関係

- ChromaticトークンのSecrets登録。
- OpenAPIが有効であること（Schemathesis）。

## ロールバック手順

- Nx/CI/Storybook関連の設定ファイルを元に戻す。
- 追加したテストを削除して従来フローに復帰。

## ヒューマンチェックポイント

- Spec承認。
- Nx構成案のレビュー。
- CI構成とChromatic連携の確認。
- リファクタリング方針のレビュー。

## 作業時間 / 担当

- 目安: 1〜2日（Nx/Storybook/CI整備が主）。
- 担当: Codex / レビュー: ユーザー。

## 実施結果

- 未実施。
