from __future__ import annotations

import pytest
from ase import Atoms as ASEAtoms

from services.supercells import build_supercell_from_grid
from app.schemas.supercell import SupercellGrid, SupercellGridAxis


def test_build_supercell_from_grid_default_axis():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    tile_a = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=base_atoms.get_cell(),
        pbc=True,
    )
    tile_b = ASEAtoms(
        symbols=["He"],
        positions=[(0.5, 0.0, 0.0)],
        cell=base_atoms.get_cell(),
        pbc=True,
    )

    grid = SupercellGrid(
        rows=2,
        cols=2,
        tiles=[["a", "b"], ["a", "b"]],
    )
    atoms_out, overlap_count = build_supercell_from_grid(
        base_atoms,
        grid,
        {"a": tile_a, "b": tile_b},
    )

    assert overlap_count == 0
    assert atoms_out.get_cell()[0][0] == 2.0
    assert atoms_out.get_cell()[1][1] == 4.0
    positions = atoms_out.get_positions()
    assert positions[0][0] == 0.0
    assert positions[1][0] == 1.5
    assert positions[2][1] == 2.0


def test_build_supercell_from_grid_axis_swap():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    grid = SupercellGrid(
        rows=3,
        cols=2,
        tiles=[["a", "a"], ["a", "a"], ["a", "a"]],
        axis=SupercellGridAxis(row="a", col="b"),
    )
    atoms_out, _ = build_supercell_from_grid(
        base_atoms,
        grid,
        {"a": base_atoms},
    )
    assert atoms_out.get_cell()[0][0] == 3.0
    assert atoms_out.get_cell()[1][1] == 4.0


def test_build_supercell_from_grid_overlap_tolerance_none_raises():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    grid = SupercellGrid(rows=1, cols=1, tiles=[["a"]])

    with pytest.raises(ValueError):
        build_supercell_from_grid(
            base_atoms,
            grid,
            {"a": base_atoms},
            check_overlap=True,
            overlap_tolerance=None,
        )


def test_build_supercell_from_grid_overlap_tolerance_non_positive_raises():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0)],
        pbc=True,
    )
    grid = SupercellGrid(rows=1, cols=1, tiles=[["a"]])

    with pytest.raises(ValueError):
        build_supercell_from_grid(
            base_atoms,
            grid,
            {"a": base_atoms},
            check_overlap=True,
            overlap_tolerance=0.0,
        )


def test_build_supercell_from_grid_overlap_count_is_per_new_atom():
    base_atoms = ASEAtoms(
        symbols=["H"],
        positions=[(0.0, 0.0, 0.0)],
        cell=[(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
        pbc=True,
    )
    tile = ASEAtoms(
        symbols=["H", "H"],
        positions=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
        cell=base_atoms.get_cell(),
        pbc=True,
    )
    grid = SupercellGrid(
        rows=1,
        cols=2,
        tiles=[["a", "a"]],
        axis=SupercellGridAxis(row="a", col="b"),
    )

    _, overlap_count = build_supercell_from_grid(
        base_atoms,
        grid,
        {"a": tile},
        check_overlap=True,
        overlap_tolerance=1.0e-6,
    )

    # 2つ目以降の各新規原子は、既存原子と1件でも重なれば+1される。
    assert overlap_count == 3
