from __future__ import annotations

from ase import Atoms as ASEAtoms

from app.schemas.supercell import SupercellGrid, SupercellGridAxis

def build_supercell_from_grid(
    base_atoms: ASEAtoms,
    grid: SupercellGrid,
    tiles: dict[str, ASEAtoms],
    *,
    check_overlap: bool = False,
    overlap_tolerance: float | None = None,
) -> tuple[ASEAtoms, int]:
    if check_overlap and (overlap_tolerance is None or overlap_tolerance <= 0):
        raise ValueError(
            "overlap_tolerance must be a positive number when check_overlap is enabled"
        )

    axis = grid.axis or SupercellGridAxis(row="b", col="a")
    base_cell = base_atoms.get_cell()
    base_a = base_cell[0]
    base_b = base_cell[1]

    row_vec = base_a if axis.row == "a" else base_b
    col_vec = base_a if axis.col == "a" else base_b

    symbols: list[str] = []
    positions: list[tuple[float, float, float]] = []
    overlap_count = 0
    tolerance_sq = (overlap_tolerance or 0.0) ** 2

    for row_index, row in enumerate(grid.tiles):
        for col_index, structure_id in enumerate(row):
            atoms = tiles[structure_id]
            shift = row_vec * row_index + col_vec * col_index
            for symbol, pos in zip(
                atoms.get_chemical_symbols(), atoms.get_positions(), strict=True
            ):
                new_pos = (
                    float(pos[0] + shift[0]),
                    float(pos[1] + shift[1]),
                    float(pos[2] + shift[2]),
                )
                if check_overlap:
                    for existing in positions:
                        dx = existing[0] - new_pos[0]
                        dy = existing[1] - new_pos[1]
                        dz = existing[2] - new_pos[2]
                        if dx * dx + dy * dy + dz * dz <= tolerance_sq:
                            overlap_count += 1
                            break
                positions.append(new_pos)
                symbols.append(symbol)

    # axis.row and axis.col are guaranteed distinct by SupercellGridAxis validation.
    repeat_a = grid.rows if axis.row == "a" else grid.cols
    repeat_b = grid.rows if axis.row == "b" else grid.cols

    out_cell = base_cell.copy()
    out_cell[0] = base_a * repeat_a
    out_cell[1] = base_b * repeat_b

    atoms_out = ASEAtoms(
        symbols=symbols,
        positions=positions,
        cell=out_cell,
        pbc=base_atoms.get_pbc(),
    )
    return atoms_out, overlap_count
