import type { Atom } from '../../lib/types'

type Vec3 = { x: number; y: number; z: number }

const filterValidIndices = (atoms: Array<Atom>, indices: Array<number>) =>
  indices.filter((index) => index >= 0 && index < atoms.length)

const calculateCentroid = (atoms: Array<Atom>) => {
  const sum = atoms.reduce(
    (acc, atom) => ({
      x: acc.x + atom.x,
      y: acc.y + atom.y,
      z: acc.z + atom.z,
    }),
    { x: 0, y: 0, z: 0 },
  )
  const inv = 1 / atoms.length
  return {
    x: sum.x * inv,
    y: sum.y * inv,
    z: sum.z * inv,
  }
}

/** 指定インデックスの原子だけを平行移動した新しい配列を返す。 */
export function shiftAtoms(
  atoms: Array<Atom>,
  indices: Array<number>,
  delta: Vec3,
): Array<Atom> {
  if (indices.length === 0) return atoms
  const set = new Set(indices)
  return atoms.map((atom, index) =>
    set.has(index)
      ? {
          ...atom,
          x: atom.x + delta.x,
          y: atom.y + delta.y,
          z: atom.z + delta.z,
        }
      : atom,
  )
}

/** 指定原子群の先頭原子を原点へ移動する。 */
export function alignSelectedToOrigin(
  atoms: Array<Atom>,
  indices: Array<number>,
): Array<Atom> {
  const validIndices = filterValidIndices(atoms, indices)
  if (validIndices.length === 0) return atoms
  const anchorIndex = validIndices[0]
  if (anchorIndex < 0 || anchorIndex >= atoms.length) return atoms
  const anchor = atoms[anchorIndex]
  return shiftAtoms(atoms, validIndices, {
    x: -anchor.x,
    y: -anchor.y,
    z: -anchor.z,
  })
}

/** 指定原子群の重心を原点へ移動する。 */
export function alignSelectedCentroid(
  atoms: Array<Atom>,
  indices: Array<number>,
): Array<Atom> {
  if (indices.length === 0) return atoms
  const selectedIndices = filterValidIndices(atoms, indices)
  if (selectedIndices.length === 0) return atoms
  const selected = selectedIndices.map((index) => atoms[index])
  const centroid = calculateCentroid(selected)
  return shiftAtoms(atoms, selectedIndices, {
    x: -centroid.x,
    y: -centroid.y,
    z: -centroid.z,
  })
}
