from __future__ import annotations

from typing import Any, Protocol, SupportsFloat, cast

from ase import Atoms as ASEAtoms

from app.schemas.supercell import SupercellGrid, SupercellGridAxis


type Vec3 = tuple[float, float, float]
type PBCLike = bool | tuple[bool, bool, bool]


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


def _extract_symbols_positions(atoms: ASEAtoms) -> list[tuple[str, Vec3]]:
    """構造の元素記号と座標を抽出する．

    Parameters:
        atoms: 抽出対象の ASE Atoms．

    Returns:
        `(symbol, position)` の配列．
    """
    symbols = cast(list[str], cast(Any, atoms).get_chemical_symbols())
    raw_positions = cast(list[object], cast(Any, atoms).get_positions())
    parsed: list[tuple[str, Vec3]] = []
    for symbol, position in zip(symbols, raw_positions, strict=True):
        parsed.append((symbol, _extract_vec3(position)))
    return parsed


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


def _add_vec3(lhs: Vec3, rhs: Vec3) -> Vec3:
    """2つの3次元ベクトルを加算する．

    Parameters:
        lhs: 左辺ベクトル．
        rhs: 右辺ベクトル．

    Returns:
        加算結果．
    """
    return (lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2])


def _scale_vec3(vec: Vec3, scalar: int) -> Vec3:
    """3次元ベクトルを整数倍する．

    Parameters:
        vec: 元のベクトル．
        scalar: 乗算係数．

    Returns:
        スカラー倍したベクトル．
    """
    return (vec[0] * scalar, vec[1] * scalar, vec[2] * scalar)


def _shifted_positions_for_tile(
    atoms: ASEAtoms,
    *,
    row_vec: Vec3,
    col_vec: Vec3,
    row_index: int,
    col_index: int,
) -> list[tuple[str, Vec3]]:
    """タイル配置に応じた平行移動後の原子座標を返す．

    Parameters:
        atoms: タイルとして配置する構造．
        row_vec: 行方向の格子ベクトル．
        col_vec: 列方向の格子ベクトル．
        row_index: 行インデックス．
        col_index: 列インデックス．

    Returns:
        平行移動後の `(symbol, position)` 配列．
    """
    shift = _add_vec3(_scale_vec3(row_vec, row_index), _scale_vec3(col_vec, col_index))
    return [
        (symbol, _add_vec3(pos, shift))
        for symbol, pos in _extract_symbols_positions(atoms)
    ]


def _count_overlap_linear(
    *, existing_positions: list[Vec3], new_pos: Vec3, tolerance_sq: float
) -> bool:
    """新規原子が既存原子と重なるかを線形探索で判定する．

    Parameters:
        existing_positions: 既存原子の座標配列．
        new_pos: 判定対象の新規原子座標．
        tolerance_sq: 許容距離の2乗値．

    Returns:
        1件以上重なりがある場合は `True`．

    Notes:
        返り値を使って `overlap_count` を1原子につき最大1回だけ増やす．
    """
    for existing in existing_positions:
        dx = existing[0] - new_pos[0]
        dy = existing[1] - new_pos[1]
        dz = existing[2] - new_pos[2]
        if dx * dx + dy * dy + dz * dz <= tolerance_sq:
            return True
    return False


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

    row_vec = base_a if axis.row == "a" else base_b
    col_vec = base_a if axis.col == "a" else base_b
    symbols: list[str] = []
    positions: list[Vec3] = []
    overlap_count = 0

    for row_index, row in enumerate(grid.tiles):
        for col_index, structure_id in enumerate(row):
            atoms = tiles[structure_id]
            shifted_atoms = _shifted_positions_for_tile(
                atoms,
                row_vec=row_vec,
                col_vec=col_vec,
                row_index=row_index,
                col_index=col_index,
            )
            for symbol, new_pos in shifted_atoms:
                if check_overlap:
                    if _count_overlap_linear(
                        existing_positions=positions,
                        new_pos=new_pos,
                        tolerance_sq=tolerance_sq,
                    ):
                        overlap_count += 1
                positions.append(new_pos)
                symbols.append(symbol)

    out_cell = _compute_output_cell(
        base_cell=base_cell,
        base_a=base_a,
        base_b=base_b,
        rows=grid.rows,
        cols=grid.cols,
        axis=axis,
    )

    atoms_out = ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=out_cell,
        pbc=_extract_pbc(base_atoms),
    )
    return atoms_out, overlap_count
