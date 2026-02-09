# Graph Visualization

報告書向けに、プロジェクト構造を次の4系統で可視化する。

- Nx graph: ワークスペース内のプロジェクト依存
- dependency-cruiser: `apps/web/src` 内のモジュール依存
- madge: `features/editor-v2/components` 内の依存関係（局所可視化）
- pyreverse: `apps/api/app` と `apps/api/services` のクラス/パッケージ依存

## Prerequisites

- 依存インストール済み: `pnpm install`
- SVG出力を行う場合のみ Graphviz の `dot` コマンドが必要

## Generate

### 1) Nx project graph

```bash
just graph-nx
```

出力:

- `docs/graphs/nx-project-graph.html`

### 2) Web module dependency graph

```bash
just graph-web-deps
```

出力:

- `docs/graphs/web-dependency-graph.dot`
- `docs/graphs/web-dependency-graph.mmd`
- `docs/graphs/web-dependency-graph.json`
- `docs/graphs/web-dependency-graph.svg` (`dot` がある場合のみ)

### 3) Editor-v2 dependency graph (madge)

```bash
just graph-editor-v2-madge
```

出力:

- `docs/graphs/editor-v2-madge.json`
- `docs/graphs/editor-v2-madge.dot`
- `docs/graphs/editor-v2-madge.circular.json`
- `docs/graphs/editor-v2-madge.svg` (`dot` がある場合のみ)

### 4) API class/package graph (pyreverse)

```bash
just api-pyreverse
```

出力:

- `docs/graphs/api-classes.dot`
- `docs/graphs/api-packages.dot`
- `docs/graphs/api-classes.svg`
- `docs/graphs/api-packages.svg`

## Notes

- 設定ファイル: `.dependency-cruiser.cjs`
- 現在は次のルールを検証する:
  - `apps/web/src/components/ui` から `apps/web/src/features` への依存を禁止
- Story / Test / build artifacts は除外対象
