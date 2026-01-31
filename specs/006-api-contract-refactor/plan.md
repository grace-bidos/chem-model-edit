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
10. Integration & Contract Testing（apps/api/tests と packages/api-client の整合テスト）。
11. Staging Deployment & Validation（API + Web をステージングへ、smoke/E2E）。
12. Performance Testing（新エンドポイントの負荷試験と閾値記録）。
13. Rollback Validation（手順のリハーサルと互換性確認）。
14. 主要フローの動作確認（parse/structures/supercell/zpe/auth）。

## 完了条件
- `/api/docs` と `/api/openapi.json` が到達可能。
- すべてのエンドポイントがOpenAPI準拠のスキーマを返す。
- TypeScriptクライアントが生成・統合されている。
- フロントに旧API型参照が残っていない。
- 主要フローのリグレッションテストが通過。
- Contract/Integrationテストが通過。

## 依存関係
- OpenAPI生成ツールの導入（openapi-typescript等）とCIジョブ追加。
- 既存クライアント/統合（apps/web、外部ZPE workerがあれば同時更新）。
- 環境別設定（API_BASE/VITE_API_BASE、CORS許可ドメイン、認証ヘッダ）。
- セキュリティ/ゲートウェイ設定（CORS/ヘッダ許可、権限トークン運用）。

## ロールバック手順
- フロントを先に旧版へ戻し、API互換性を確保。
- `main.py` と旧ルーティングを復元し、新APIを無効化。
- 生成クライアント/型参照を旧版へ差し戻す。
- （DB変更が入った場合）マイグレーションのロールバック or ガード適用。
- 監視/アラートを確認し、必要に応じてキャッシュ/状態をクリア。

## ヒューマンチェックポイント
- Spec承認。
- API URL/契約レビュー。
- 生成クライアントのAPI表現レビュー。

## 作業時間 / 担当
- 目安: 1〜2週間。
- Phase 1: Backend契約定義とOpenAPI整備（2〜3日）。
- Phase 2: クライアント生成とフロント統合（2〜3日）。
- Phase 3: テスト/検証/調整（2〜3日）。
- 担当: Codex / レビュー: ユーザー。

## 実施結果
- 未実施。
