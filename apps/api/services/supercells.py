from __future__ import annotations

from typing import Any, Protocol, SupportsFloat, cast

import numpy as np
from numpy.typing import NDArray
from ase import Atoms as ASEAtoms

from app.schemas.supercell import SupercellGrid, SupercellGridAxis


type Vec3 = tuple[float, float, float]
type PBCLike = bool | tuple[bool, bool, bool]
type FloatArray = NDArray[np.float64]


class CellLike(Protocol):
    def copy(self) -> "CellLike":
        ...

    def __setitem__(self, key: int, value: Vec3) -> None:
        ...


def _validate_overlap_options(
    *, check_overlap: bool, overlap_tolerance: float | None
) -> float:
    """重なり判定オプションを検証し，2乗許容距離を返す．

    Parameters:
        check_overlap: 重なり判定を有効化するかどうか．
        overlap_tolerance: 重なり判定の許容距離．

    Returns:
        重なり判定に使う許容距離の2乗値．

    Raises:
        ValueError: `check_overlap=True` かつ `overlap_tolerance` が不正な場合．
    """
    if check_overlap and (overlap_tolerance is None or overlap_tolerance <= 0):
        raise ValueError(
            "overlap_tolerance must be a positive number when check_overlap is enabled"
        )
    return (overlap_tolerance or 0.0) ** 2


def _resolve_axis(axis: SupercellGridAxis | None) -> SupercellGridAxis:
    """グリッド軸マッピングを解決する．

    Parameters:
        axis: リクエストで指定された軸マッピング．

    Returns:
        解決済みの軸マッピング．未指定時は `row=b, col=a`．
    """
    return axis or SupercellGridAxis(row="b", col="a")


def _extract_vec3(value: object) -> Vec3:
    """配列ライクな3次元ベクトルを `Vec3` へ正規化する．

    Parameters:
        value: 3要素の数値配列．

    Returns:
        `Vec3` に変換済みの座標．
    """
    vec = cast(tuple[object, object, object], value)
    return (
        float(cast(SupportsFloat, vec[0])),
        float(cast(SupportsFloat, vec[1])),
        float(cast(SupportsFloat, vec[2])),
    )


def _extract_cell_vectors(base_atoms: ASEAtoms) -> tuple[Vec3, Vec3, CellLike]:
    """基準構造から格子ベクトルを取り出す．

    Parameters:
        base_atoms: 基準となる ASE Atoms．

    Returns:
        `a` ベクトル，`b` ベクトル，出力セル作成に使うセルオブジェクト．
    """
    base_cell = cast(CellLike, cast(Any, base_atoms).get_cell())
    base_a = _extract_vec3(cast(Any, base_cell)[0])
    base_b = _extract_vec3(cast(Any, base_cell)[1])
    return base_a, base_b, base_cell


def _extract_pbc(base_atoms: ASEAtoms) -> PBCLike:
    """周期境界条件を抽出する．

    Parameters:
        base_atoms: 基準となる ASE Atoms．

    Returns:
        ASE Atoms 生成に渡せる周期境界条件．
    """
    raw_pbc = cast(Any, base_atoms).get_pbc()
    if isinstance(raw_pbc, bool):
        return raw_pbc
    pbc = cast(tuple[object, object, object], raw_pbc)
    return (bool(pbc[0]), bool(pbc[1]), bool(pbc[2]))


def _extract_symbols_positions_array(atoms: ASEAtoms) -> tuple[list[str], FloatArray]:
    """構造から元素記号配列と座標配列を抽出する．

    Parameters:
        atoms: 抽出対象の ASE Atoms．

    Returns:
        元素記号配列と，`(N, 3)` 形状の座標配列．
    """
    symbols = cast(list[str], cast(Any, atoms).get_chemical_symbols())
    positions = cast(Any, atoms).get_positions()
    pos_array = np.asarray(positions, dtype=np.float64)
    return symbols, pos_array


def _scale_vec3(vec: Vec3, scalar: int) -> Vec3:
    """3次元ベクトルを整数倍する．

    Parameters:
        vec: 元のベクトル．
        scalar: 乗算係数．

    Returns:
        スカラー倍したベクトル．
    """
    return (vec[0] * scalar, vec[1] * scalar, vec[2] * scalar)


def _compute_output_cell(
    *,
    base_cell: CellLike,
    base_a: Vec3,
    base_b: Vec3,
    rows: int,
    cols: int,
    axis: SupercellGridAxis,
) -> CellLike:
    """出力スーパーセルの格子を計算する．

    Parameters:
        base_cell: 基準構造のセル．
        base_a: 基準構造の a ベクトル．
        base_b: 基準構造の b ベクトル．
        rows: グリッド行数．
        cols: グリッド列数．
        axis: 行列と格子軸の対応．

    Returns:
        拡張済みのセルオブジェクト．
    """
    # axis.row and axis.col are guaranteed distinct by SupercellGridAxis validation.
    repeat_a = rows if axis.row == "a" else cols
    repeat_b = rows if axis.row == "b" else cols

    out_cell = base_cell.copy()
    out_cell[0] = _scale_vec3(base_a, repeat_a)
    out_cell[1] = _scale_vec3(base_b, repeat_b)
    return out_cell


def _count_incremental_overlaps_np(
    *,
    existing_positions: FloatArray,
    new_positions: FloatArray,
    tolerance_sq: float,
) -> int:
    """新規配置原子に対して，既存原子との重なり数を増分計算する．

    Parameters:
        existing_positions: すでに配置済みの座標配列．
        new_positions: 今回追加する座標配列．
        tolerance_sq: 許容距離の2乗値．

    Returns:
        新規原子に対する重なりヒット数．
    """
    overlap_count = 0
    for index in range(new_positions.shape[0]):
        point = new_positions[index]
        hit_existing = False
        if existing_positions.shape[0] > 0:
            diff_existing = existing_positions - point
            d2_existing = np.einsum("ij,ij->i", diff_existing, diff_existing)
            hit_existing = bool(np.any(d2_existing <= tolerance_sq))
        if hit_existing:
            overlap_count += 1
            continue
        if index > 0:
            prev = new_positions[:index]
            diff_prev = prev - point
            d2_prev = np.einsum("ij,ij->i", diff_prev, diff_prev)
            if bool(np.any(d2_prev <= tolerance_sq)):
                overlap_count += 1
    return overlap_count


def count_overlap_pairs_spatial_hash(points: FloatArray, *, tolerance: float) -> int:
    """空間ハッシュ法で閾値内ペア数を数える．

    Parameters:
        points: `(N, 3)` 形状の座標配列．
        tolerance: 判定距離．

    Returns:
        閾値内にある原子ペア数．

    Raises:
        ValueError: `tolerance <= 0` の場合．
    """
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    if points.shape[0] < 2:
        return 0

    tolerance_sq = tolerance * tolerance
    cell_size = tolerance
    buckets: dict[tuple[int, int, int], list[int]] = {}
    count = 0

    offsets = [
        (dx, dy, dz)
        for dx in (-1, 0, 1)
        for dy in (-1, 0, 1)
        for dz in (-1, 0, 1)
    ]

    for idx in range(points.shape[0]):
        point = points[idx]
        cell = tuple(np.floor(point / cell_size).astype(np.int64).tolist())
        cx, cy, cz = cell
        for ox, oy, oz in offsets:
            neighbor = (cx + ox, cy + oy, cz + oz)
            for other_idx in buckets.get(neighbor, []):
                diff = points[other_idx] - point
                d2 = float(np.dot(diff, diff))
                if d2 <= tolerance_sq:
                    count += 1
        buckets.setdefault(cell, []).append(idx)
    return count


def build_supercell_from_grid(
    base_atoms: ASEAtoms,
    grid: SupercellGrid,
    tiles: dict[str, ASEAtoms],
    *,
    check_overlap: bool = False,
    overlap_tolerance: float | None = None,
) -> tuple[ASEAtoms, int]:
    """タイルグリッド定義からスーパーセルを構築する．

    Parameters:
        base_atoms: 平行移動の基準となる格子を持つ構造．
        grid: グリッド定義．`tiles` の配置と軸対応を含む．
        tiles: `structure_id -> ASEAtoms` の参照表．
        check_overlap: 重なり判定を有効化するかどうか．
        overlap_tolerance: 重なり判定の許容距離．

    Returns:
        構築したスーパーセルと重なりカウントのタプル．

    Raises:
        ValueError: `check_overlap=True` かつ `overlap_tolerance` が不正な場合．
        KeyError: `tiles` に存在しない `structure_id` が指定された場合．
    """
    tolerance_sq = _validate_overlap_options(
        check_overlap=check_overlap,
        overlap_tolerance=overlap_tolerance,
    )
    axis = _resolve_axis(grid.axis)
    base_a, base_b, base_cell = _extract_cell_vectors(base_atoms)

    row_vec = np.asarray(base_a if axis.row == "a" else base_b, dtype=np.float64)
    col_vec = np.asarray(base_a if axis.col == "a" else base_b, dtype=np.float64)

    symbols: list[str] = []
    position_chunks: list[FloatArray] = []
    overlap_count = 0
    existing_positions = np.empty((0, 3), dtype=np.float64)

    tile_cache: dict[str, tuple[list[str], FloatArray]] = {}

    for row_index, row in enumerate(grid.tiles):
        for col_index, structure_id in enumerate(row):
            atoms = tiles[structure_id]
            cached = tile_cache.get(structure_id)
            if cached is None:
                cached = _extract_symbols_positions_array(atoms)
                tile_cache[structure_id] = cached
            tile_symbols, tile_positions = cached

            shift = row_vec * row_index + col_vec * col_index
            shifted_positions = tile_positions + shift

            if check_overlap:
                overlap_count += _count_incremental_overlaps_np(
                    existing_positions=existing_positions,
                    new_positions=shifted_positions,
                    tolerance_sq=tolerance_sq,
                )

            symbols.extend(tile_symbols)
            position_chunks.append(shifted_positions)
            if existing_positions.shape[0] == 0:
                existing_positions = shifted_positions.copy()
            else:
                existing_positions = cast(
                    FloatArray,
                    np.concatenate((existing_positions, shifted_positions), axis=0),
                )

    out_cell = _compute_output_cell(
        base_cell=base_cell,
        base_a=base_a,
        base_b=base_b,
        rows=grid.rows,
        cols=grid.cols,
        axis=axis,
    )

    if position_chunks:
        positions_array = cast(FloatArray, np.concatenate(position_chunks, axis=0))
    else:
        positions_array = np.empty((0, 3), dtype=np.float64)

    atoms_out = ASEAtoms(
        symbols=symbols,
        positions=positions_array.tolist(),
        cell=out_cell,
        pbc=_extract_pbc(base_atoms),
    )
    return atoms_out, overlap_count
