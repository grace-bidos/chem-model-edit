# Plan: API契約/URL整理とフロント連携リファクタ

## 目的
API契約を固定し、フロント/バック間の責務分割を明確化する。

## 作業ステップ
1. Specレビュー（ヒューマン承認）。
2. FastAPI構成の再編（app/core, app/schemas, app/routers, app/services, app/deps）。
3. PydanticスキーマをcamelCase対応で整理（alias/serialization設定）。
4. `/api` 配下に新ルーティングを実装（旧URLは削除）。
5. 例外ハンドラを整備し `ErrorResponse` に統一。
6. OpenAPI生成フローを追加し `packages/api-client` を構築。
7. フロントのAPIアクセスを生成クライアントへ移行。
8. 既存の手書きAPI型を削除/移動し参照を整理。
9. APIテスト/フロント型参照の更新（apps/api/tests, apps/web）。
10. 主要フローの動作確認（parse/structures/supercell/zpe/auth）。

## 完了条件
- Specの承認。
- 受け入れ基準を満たす。

## 依存関係
- OpenAPI生成ツールの導入（openapi-typescript等）。
- ZPE workerが外部にある場合は同時更新が必要。

## ロールバック手順
- `main.py` と旧ルーティングを復元。
- 生成クライアントの参照を元に戻す。
- `apps/web` のAPI呼び出しを旧版に差し戻す。

## ヒューマンチェックポイント
- Spec承認。
- API URL/契約レビュー。
- 生成クライアントのAPI表現レビュー。

## 作業時間 / 担当
- 目安: 1〜2日。
- 担当: Codex / レビュー: ユーザー。

## 実施結果
- 未実施。
