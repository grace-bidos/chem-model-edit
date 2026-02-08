import { describe, expect, it } from 'vitest'

import {
  alignSelectedCentroid,
  alignSelectedToOrigin,
  shiftAtoms,
} from './align'

import type { Atom } from '@/lib/types'

const atoms: Array<Atom> = [
  { symbol: 'H', x: 1, y: 0, z: 0 },
  { symbol: 'O', x: 3, y: 0, z: 0 },
  { symbol: 'H', x: 5, y: 0, z: 0 },
]

describe('align utilities', () => {
  it('shifts only selected atoms', () => {
    const shifted = shiftAtoms(atoms, [0, 2], { x: 1, y: -1, z: 2 })

    expect(shifted[0]).toMatchObject({ x: 2, y: -1, z: 2 })
    expect(shifted[1]).toMatchObject({ x: 3, y: 0, z: 0 })
    expect(shifted[2]).toMatchObject({ x: 6, y: -1, z: 2 })
  })

  it('aligns selected atoms to origin using first selected atom', () => {
    const aligned = alignSelectedToOrigin(atoms, [1, 2])

    expect(aligned[1]).toMatchObject({ x: 0, y: 0, z: 0 })
    expect(aligned[2]).toMatchObject({ x: 2, y: 0, z: 0 })
  })

  it('aligns selected atoms by centroid using only valid indices', () => {
    const aligned = alignSelectedCentroid(atoms, [0, 1, 10])

    expect(aligned[0].x).toBeCloseTo(-1)
    expect(aligned[1].x).toBeCloseTo(1)
    expect(aligned[2].x).toBeCloseTo(5)
  })
})
