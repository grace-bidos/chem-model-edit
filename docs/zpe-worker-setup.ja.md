# 計算ノードセットアップ（runtime専用）

この経路は Redis 非依存です。  
計算ノード登録は `/api/runtime/*` と Convex 管理のノード情報で行います。

## 前提

- control-plane API（`/api/runtime/nodes/*`）へ到達できる
- Web SaaS 側で Clerk サインインできる
- Ubuntu計算ノードに `curl` と `jq` がある

## 運用フロー

1. Web SaaS でサインインし、計算ノード追加コマンドを生成する。
2. 表示されたコマンドをコピーする。
3. 計算ノードのCLIで実行する。

生成されるコマンド形式:

```bash
curl -fsSL <API_BASE>/api/runtime/nodes/install.sh | bash -s -- \
  --api-base <API_BASE> \
  --join-token <JOIN_TOKEN> \
  --queue-name <QUEUE_NAME>
```

インストーラは `POST /api/runtime/nodes/register` を呼び出し、  
コマンドを発行したユーザーの runtime target を作成します。

## 補足

- join token は短寿命・単回利用です。
- 有効ターゲット選択は `/api/runtime/targets` で管理します。
