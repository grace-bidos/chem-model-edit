import { describe, expect, it } from 'vitest'

import {
  alignSelectedCentroid,
  alignSelectedToOrigin,
  shiftAtoms,
} from './align'

import type { Atom } from '@/lib/types'

const atoms: Array<Atom> = [
  { symbol: 'H', x: 1, y: 2, z: 3 },
  { symbol: 'O', x: 3, y: 4, z: 5 },
  { symbol: 'H', x: 5, y: 6, z: 7 },
]

describe('align utilities', () => {
  it('shifts only selected atoms', () => {
    const shifted = shiftAtoms(atoms, [0, 2], { x: 1, y: -1, z: 2 })

    expect(shifted[0]).toMatchObject({ x: 2, y: 1, z: 5 })
    expect(shifted[1]).toMatchObject({ x: 3, y: 4, z: 5 })
    expect(shifted[2]).toMatchObject({ x: 6, y: 5, z: 9 })
  })

  it('aligns selected atoms to origin using first selected atom', () => {
    const aligned = alignSelectedToOrigin(atoms, [1, 2])

    expect(aligned[1]).toMatchObject({ x: 0, y: 0, z: 0 })
    expect(aligned[2]).toMatchObject({ x: 2, y: 2, z: 2 })
  })

  it('keeps index 0 as a valid anchor when aligning to origin', () => {
    const aligned = alignSelectedToOrigin(atoms, [0, 2])

    expect(aligned[0]).toMatchObject({ x: 0, y: 0, z: 0 })
    expect(aligned[2]).toMatchObject({ x: 4, y: 4, z: 4 })
  })

  it('aligns selected atoms by centroid using only valid indices', () => {
    const aligned = alignSelectedCentroid(atoms, [0, 1, 10])

    expect(aligned[0]).toMatchObject({ x: -1, y: -1, z: -1 })
    expect(aligned[1]).toMatchObject({ x: 1, y: 1, z: 1 })
    expect(aligned[2]).toMatchObject({ x: 5, y: 6, z: 7 })
  })

  it('filters out negative and out-of-range indices for centroid alignment', () => {
    const aligned = alignSelectedCentroid(atoms, [-1, 2, 99])

    expect(aligned[2]).toMatchObject({ x: 0, y: 0, z: 0 })
    expect(aligned[0]).toMatchObject({ x: 1, y: 2, z: 3 })
    expect(aligned[1]).toMatchObject({ x: 3, y: 4, z: 5 })
  })

  it('returns original atoms when selected indices are empty or invalid', () => {
    expect(alignSelectedToOrigin(atoms, [])).toBe(atoms)
    expect(alignSelectedToOrigin(atoms, [999])).toBe(atoms)
    expect(alignSelectedCentroid(atoms, [])).toBe(atoms)
    expect(alignSelectedCentroid(atoms, [999])).toBe(atoms)
  })
})
