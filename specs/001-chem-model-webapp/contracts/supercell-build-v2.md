# Contract: Supercell Build v2 (Editor v2)

## Endpoints

- POST /structures
- POST /structures/import
- POST /supercell/build

## POST /structures (register Structure)

Request:
```json
{
  "structure": {
    "atoms": [{ "symbol": "O", "x": 0, "y": 0, "z": 0 }],
    "lattice": {
      "a": { "x": 5, "y": 0, "z": 0 },
      "b": { "x": 0, "y": 5, "z": 0 },
      "c": { "x": 0, "y": 0, "z": 5 }
    }
  },
  "source": "editor-v2"
}
```

Response:
```json
{
  "structureId": "9f1c2d7d0c9a4a6c8c5c34f3b8b5b6e3",
  "source": "editor-v2"
}
```

## POST /structures/import (QE import)

Request:
```json
{
  "content": "<QE .in>"
}
```

Response:
```json
{
  "structureId": "a6b5c4d3e2f1",
  "structure": {
    "atoms": [{ "symbol": "O", "x": 0, "y": 0, "z": 0 }],
    "lattice": {
      "a": { "x": 5, "y": 0, "z": 0 },
      "b": { "x": 0, "y": 5, "z": 0 },
      "c": { "x": 0, "y": 0, "z": 5 }
    }
  },
  "source": "ase"
}
```

## POST /supercell/build

Request:
```json
{
  "baseStructureId": "s-1",
  "grid": {
    "rows": 2,
    "cols": 3,
    "tiles": [
      ["s-1", "s-2", "s-1"],
      ["s-3", "s-1", "s-2"]
    ],
    "axis": { "row": "b", "col": "a" }
  },
  "options": {
    "checkOverlap": true,
    "overlapTolerance": 0.2,
    "validateLattice": "none"
  },
  "output": {
    "register": true,
    "includeStructure": false
  }
}
```

Response:
```json
{
  "structureId": "s-out-1",
  "meta": {
    "rows": 2,
    "cols": 3,
    "tileCount": 6,
    "overlapCount": 1,
    "baseStructureId": "s-1",
    "structureIdsUsed": ["s-1", "s-2", "s-3"]
  }
}
```

## Notes

- Unit: angstrom (Cartesian coordinates and lattice vectors).
- baseStructureId の lattice を基準にタイルを平行移動する。
- grid.axis が省略された場合は row->b, col->a をデフォルトとする。
- tiles は矩形であること、rows/cols と一致すること。
- structureId はバックエンドの StructureRegistry（ASE Atoms）を参照する。
- output.includeStructure=true の場合は Response に Structure を含める。
