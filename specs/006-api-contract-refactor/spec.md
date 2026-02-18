# Feature Specification: API契約/URL整理とフロント連携リファクタ

> [!NOTE]
> Historical spec document. Superseded by the Modal-first runtime/deploy operating model on February 18, 2026.

**Feature Branch**: `[refactor/api-contract]`  
**Created**: 2026-01-31  
**Status**: Draft  
**Input**: ユーザー要望「APIスキーマ/エンドポイントを整理し、フロント・バックの責務分割を容易にする」

## 背景 / 目的

- 現状は `apps/api/main.py` にルートが集中し、スキーマ/URLがばらけている。
- フロント側は `apps/web/src/lib/types.ts` と `packages/shared` に型が二重化し、API変更が波及しやすい。
- 今後の中身のリファクタを進めるために、**契約（API）を固定**し、フロント/バックの責務分割を明確にしたい。

## スコープ

- `/api` プレフィックスでAPIを統一（バージョン無し）。
- URL命名・構造の整理（名詞/複数形、リソース指向）。
- すべてのAPIスキーマをPydanticに集約し、OpenAPIを正とする。
- OpenAPIからTS型/クライアントを生成し、フロントは生成物のみ参照。
- フロント側のAPI呼び出しと型参照を更新。
- ルーティングを機能別に分割（routers/services/schemasなど）。

## 非ゴール

- 物理計算ロジック/アルゴリズムの変更。
- 計算結果やUI動作の変更（最終的な挙動は維持）。
- 永続ストレージの導入やDB移行。
- 外部連携方式の刷新（Cloud Run / Cloudflare Pagesの変更はしない）。

## 方針 / 決定事項

- **Base Path**: `/api` に統一（`/api/v1` は使わない）。
- **命名規則**: ワイヤ上は camelCase。内部は snake_case 可（Pydanticのaliasで変換）。
- **認証**: `Authorization: Bearer <token>` を統一。
- **エラー形式**: `ErrorResponse { error: { code, message, details? } }` を共通化。
- **OpenAPI**: `/api/openapi.json`、Docs `/api/docs` を有効化。
- **ID/時刻**: `id` は文字列、時刻は ISO 8601。
- **レスポンスID**: 新規作成やリソース応答は `id` を基本とし、クロス参照は `structureId` 等を使う。

## APIエンドポイント（新）

### 共通

- `GET /api/health`

### Auth Schemas

- `POST /api/auth/register`
- `POST /api/auth/login`
- `POST /api/auth/logout`
- `GET /api/auth/me`

### Structures Schemas

- `POST /api/structures/parse`
- `POST /api/structures`
- `GET /api/structures/{id}`
- `GET /api/structures/{id}/view?format=cif`
- `POST /api/structures/export`

### Transforms / Lattice Schemas

- `POST /api/transforms/delta-transplant`
- `POST /api/lattices/convert`

### Supercell Schemas

- `POST /api/supercells`
- `POST /api/supercells/tiled`
- `POST /api/supercells/builds`

### ZPE (user) Schemas

- `POST /api/zpe/parse`
- `POST /api/zpe/jobs`
- `GET /api/zpe/jobs/{id}`
- `GET /api/zpe/jobs/{id}/result`
- `GET /api/zpe/jobs/{id}/files?kind=summary|freqs`

### ZPE (compute/admin) Schemas

- `POST /api/zpe/compute/enroll-tokens`
- `POST /api/zpe/compute/servers`
- `DELETE /api/zpe/compute/servers/{id}`
- `POST /api/zpe/compute/jobs/lease`
- `POST /api/zpe/compute/jobs/{id}/result`
- `POST /api/zpe/compute/jobs/{id}/failed`
- `GET /api/zpe/targets?limit=50&offset=0`
- `PUT /api/zpe/targets/{id}/active`
- `GET /api/zpe/admin/ops`
- `PATCH /api/zpe/admin/ops`

## 契約（スキーマ）

### 共通型

- `Atom`: `{ symbol, x, y, z }`
- `Vector3`: `{ x, y, z }`
- `Lattice`: `{ a, b, c }`
- `LatticeParams`: `{ a, b, c, alpha, beta, gamma }`
- `Structure`: `{ atoms, lattice? }`
- `QeParameters`: `{ control, system, electrons, ions, cell, pseudopotentials, kpoints? }`
- `Pagination`: `{ total, limit, offset }`

### Auth

- `AuthRegisterRequest`: `{ email, password }`
- `AuthLoginRequest`: `{ email, password }`
- `AuthUser`: `{ id, email, createdAt }`
- `AuthSession`: `{ token, expiresAt, user }`
- `AuthMe`: `{ user, expiresAt }`
- `AuthLogoutResponse`: `{ ok: true }`

### Structures

- `StructureParseRequest`: `{ content, format? }`
- `StructureParseResponse`: `{ structure }`
- `StructureCreateRequest`: `{ content, format? }`
- `StructureCreateResponse`: `{ id, structure, source, params?, rawInput? }`
- `StructureGetResponse`: `{ structure, params?, rawInput?, source? }`
- `StructureExportRequest`: `{ structure, format? }`
- `StructureExportResponse`: `{ content }`
  - `/api/structures/{id}/view` は `chemical/x-cif` を返す（JSONではない）。

### Transforms / Lattice

- `DeltaTransplantRequest`: `{ smallIn, smallOut, largeIn }`
- `DeltaTransplantResponse`: `{ content }`
- `LatticeConvertRequest`: `{ from, unit?, lattice?, params? }`
  - `from` は `"vectors" | "params"`、対応する入力フィールド必須。
  - `unit` の既定値は `"angstrom"`。
- `LatticeConvertResponse`: `{ lattice, params, unit }`

### Supercell

- `SupercellRequest`: `{ structureA, structureB, sequence, lattice }`
- `SupercellResponse`: `{ structure, meta }`
- `SupercellBuildRequest`: `{ baseStructureId, grid, options?, output? }`
- `SupercellBuildResponse`: `{ id, structure?, meta }`

### ZPE (user)

- `ZPEParseRequest`: `{ content, structureId? }`
- `ZPEParseResponse`: `{ structure, fixedIndices, atomicSpecies, kpoints? }`
- `ZPEJobRequest`: `{ content, mobileIndices, useEnviron?, inputDir?, calcMode?, structureId? }`
  - `calcMode` は `"new" | "continue"`、既定値は `"continue"`。
- `ZPEJobResponse`: `{ id }`
- `ZPEJobStatus`: `{ status, detail?, updatedAt? }`
- `ZPEJobResultResponse`: `{ result }`
  - `/api/zpe/jobs/{id}/files` は `text/plain`（summary）または `text/csv`（freqs）を返す。
- `ZPEResult`: `{ freqsCm, zpeEv, sVibJmolK, mobileIndices, fixedIndices, kpoints, delta, lowCutCm, temperature, useEnviron, qeInput, pseudoDir, calcStartTime, calcEndTime, elapsedSeconds, cacheChecked, cacheDeleted, ecutwfc?, ecutrho? }`

### ZPE (compute/admin)

- `EnrollTokenRequest`: `{ ttlSeconds?, label? }`
- `EnrollTokenResponse`: `{ token, expiresAt, ttlSeconds, label? }`
- `ComputeRegisterRequest`: `{ token, name?, queueName?, meta? }`
- `ComputeRegisterResponse`: `{ id, registeredAt, name?, workerToken, tokenExpiresAt, tokenTtlSeconds }` (id は serverId)
- `ComputeRevokeResponse`: `{ revokedCount }`
- `QueueTarget`: `{ id, queueName, serverId, registeredAt, name? }`
- `QueueTargetListResponse`: `{ targets, activeTargetId?, pagination }`
  - `pagination`: `{ total, limit, offset }`（既定: `limit=50`, `offset=0`）
- `QueueTargetSelectResponse`: `{ activeTargetId }`
- `ComputeLeaseResponse`: `{ jobId, payload, leaseId, leaseTtlSeconds, meta? }`
- `ComputeResultRequest`: `{ leaseId, result, summaryText, freqsCsv, meta? }`
- `ComputeResultResponse`: `{ ok, idempotent }`
- `ComputeFailedRequest`: `{ leaseId, errorCode, errorMessage, traceback? }`
- `ComputeFailedResponse`: `{ ok, requeued, retryCount }`
- `OpsFlagsRequest`: `{ submissionEnabled?, dequeueEnabled? }`
- `OpsFlagsResponse`: `{ submissionEnabled, dequeueEnabled }`
  - `GET /api/zpe/admin/ops` は単一のフラグオブジェクトを返す（リストではない）。

### エラー

- `ErrorResponse`: `{ error: { code, message, details? } }`
  - すべての 4xx/5xx で統一形式にする（FastAPIの例外ハンドラを差し替え）。
  - `code` は `not_found` / `validation_error` / `unauthorized` / `conflict` などの機械可読値。

## フロントエンド変更

- OpenAPI生成クライアントを `packages/api-client` に追加。
- `apps/web` は `@chem-model/api-client` の型とリクエスト関数のみ利用。
- `apps/web/src/lib/types.ts` のAPI型は削除 or 移動。
- `apps/web/src/lib/api.ts` を生成クライアントに置き換え。
- `apps/web/src/server/api.ts` のエラーハンドリングを `ErrorResponse` に合わせる。
- `API_BASE` は `/api` 前提で解決（`http://localhost:8000/api` など）。
- 既存UIの動作は維持（見た目/操作は変えない）。

## 成果物（期待される状態）

1. `/api` 配下に統一されたエンドポイント。
2. Pydanticスキーマが単一の契約源となる。
3. OpenAPIから生成されたTS型/クライアントが`packages/api-client`に存在。
4. フロントは生成物に依存し、手書きAPI型を持たない。
5. ルーティングが `routers/` に分割され、責務が明確になる。

## 受け入れ基準

- `/api/openapi.json` と `/api/docs` が確認できる。
- 主要フロー（parse, structures, supercell, zpe）が従来と同じ動作。
- フロントのAPI型が `@chem-model/api-client` に統一される。
- 旧エンドポイントへの依存がなくなる。

## 移行方針

- 旧APIは即時廃止（互換維持は行わない）。
- フロント/バックは同時更新を前提とする（外部クライアントが無い前提）。
- 必要なら旧URLの簡易リダイレクトは検討（最小限）。

## リスク / 未確定事項

- 外部のZPE workerが存在する場合、同時更新が必要。
- camelCase化により既存のフロント/バッチが壊れる可能性。
- 生成クライアント導入に伴うビルド/CI設定の微調整。
