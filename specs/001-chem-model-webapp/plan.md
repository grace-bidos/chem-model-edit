# Implementation Plan: 計算化学向け構造編集Webアプリ

**Branch**: `[001-chem-model-webapp]` | **Date**: 2026-01-05 | **Spec**: `specs/001-chem-model-webapp/spec.md`
**Input**: Feature specification from `specs/001-chem-model-webapp/spec.md`

## Summary

Quantum ESPRESSO の .in を中心に、2構造(A/B)の可視化・編集・表面転写・複合スーパーセル構築・距離QC・格子変換を行うWebアプリを構築する。フロントは TanStack Start + shadcn/ui + Mol*、バックエンドは FastAPI + ASE/pymatgen に分離する。P1（読み込み/編集/表示）を最優先で実装し、P2（表面転写/距離QC）とP3（タイル合成/格子変換）を段階的に追加する。

## Technical Context

**Language/Version**: TypeScript (React), Python 3.11+  
**Primary Dependencies**: TanStack Start, shadcn/ui, Mol*, FastAPI, ASE, pymatgen  
**Storage**: 基本はファイル/メモリ  
**Testing**: Frontend: Vitest + Playwright(予定), Backend: pytest  
**Target Platform**: モダンブラウザ（Chrome/Firefox/Safari最新系） + ローカル/サーバ運用  
**Project Type**: Web (frontend + backend)  
**Performance Goals**: ~1k原子で 60fps 近い操作感、編集反映 <300ms  
**Constraints**: 対象は QE .in、出力は入力形式/単位を保持がデフォルト、celldmは保持（再計算/削除は将来対応）  
**Scale/Scope**: 2画面（Editor / Supercell）中心の中規模、A/Bの2構造を前提

## Constitution Check

- `.specify/memory/constitution.md` が未整備のためゲートは未定義。必要なら先に合意済み原則を記載する。

## Project Structure

### Documentation (this feature)

```text
specs/001-chem-model-webapp/
├── spec.md
├── plan.md
├── research.md        # (必要に応じて追加)
├── data-model.md      # (必要に応じて追加)
├── quickstart.md      # (必要に応じて追加)
├── contracts/         # (API契約を追加)
└── tasks.md           # (Taskステージで作成)
```

### Source Code (repository root)

```text
apps/
├── web/               # TanStack Start + shadcn/ui + Mol*
└── api/               # FastAPI + ASE/pymatgen

packages/
└── shared/            # 型定義/ユーティリティ（必要なら）
```

**Structure Decision**: フロント/バック分離を前提に monorepo 構成。初期は apps/web と apps/api を用意し、共通型を packages/shared に配置する計画。

## Phased Plan

### Phase 0: 仕様確定（Spec承認）
- Spec のレビューと修正
- 出力オプションは設定画面で管理し、デフォルトは入力形式/単位保持を確定
- celldm は保持のみ（再計算/削除は将来対応）を明記
- 取り扱う QE .in の範囲は当面「原子種・座標・格子」に限定する

### Phase 1: 基盤設計
- データモデル設計（Model(A/B)/Structure/Atom/Lattice/Selection/SurfaceTransfer/DistanceReport/TilingPattern/ExportSettings）
- FastAPI API契約（parse, export, transfer, distance, supercell, lattice-convert）
- Mol* 表示統合設計（入力フォーマット変換、表示/色/ラベル/背景）

### Phase 2: MVP実装 (P1)
- .in 読み込み/パース（ASE/pymatgen）
- 3D表示 + テーブル編集の同期
- 基本編集操作（座標/元素変更、コピー/貼り付け）
- エクスポート（QE .in への書き戻し、入力形式/単位保持）
- 出力設定画面（デフォルト保持のみ、将来拡張を見据えたUI枠）

### Phase 3: 表面転写/距離QC (P2)
- 2構造(A/B)の並列表示
- 選択UI（表面領域の選択を含む）
- インデックス対応の表面転写（Cartesian座標上書き）
- 対応原子距離の一覧（PBC時は最小像）
- A/B不整合時の事前警告（原子数・格子非互換）

### Phase 4: 複合スーパーセル/格子変換 (P3)
- タイルパターン入力（2D配列 or GUI）
- タイル合成（タイル内座標保持＋並進のみ）
- 重複/衝突チェック（許容誤差つき、任意）
- 格子パラメータ変換（格子定数 ⇄ セルベクトル）

## Risks & Mitigations

- **QE .in の多様性**: ASE/pymatgen の対応範囲を明確化し、非対応は明示的エラーにする
- **Mol* への変換**: PDB/CIF 中間形式を採用し、変換処理を共通化
- **インデックス対応の前提**: A/Bの整列・並び順が崩れると誤転写の恐れがあるため、事前警告とガイドを用意
- **celldm 整合性**: 出力保持の前提で、セル編集時の不整合に注意喚起を出す

## Open Questions

- タイルパターン入力のUI設計（2D配列/GUIの詳細）
- 表面選択の補助（Z閾値など）を導入するか
- celldm の再計算/削除をどの段階で実装するか
- 取り扱う QE .in の詳細仕様（constraints, magnetization, tags など）は後続
