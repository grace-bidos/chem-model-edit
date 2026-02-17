import fc from 'fast-check'
import { describe, expect, it } from 'vitest'

import { shiftAtoms } from './align'

import type { Atom } from '@/lib/types'

const atomArb = fc.record({
  symbol: fc.stringMatching(/^[A-Z][a-z]{0,2}$/),
  x: fc.float({ min: -1000, max: 1000, noNaN: true }),
  y: fc.float({ min: -1000, max: 1000, noNaN: true }),
  z: fc.float({ min: -1000, max: 1000, noNaN: true }),
}) satisfies fc.Arbitrary<Atom>

describe('shiftAtoms properties', () => {
  it('preserves non-selected atoms and shifts selected atoms by delta', () => {
    fc.assert(
      fc.property(
        fc.array(atomArb, { minLength: 1, maxLength: 30 }),
        fc.array(fc.nat({ max: 29 }), { minLength: 1, maxLength: 30 }),
        fc.record({
          x: fc.float({ min: -10, max: 10, noNaN: true }),
          y: fc.float({ min: -10, max: 10, noNaN: true }),
          z: fc.float({ min: -10, max: 10, noNaN: true }),
        }),
        (atoms, indices, delta) => {
          const bounded = indices
            .map((index) => index % atoms.length)
            .filter((index, i, arr) => arr.indexOf(index) === i)
          const selected = new Set(bounded)
          const shifted = shiftAtoms(atoms, bounded, delta)

          for (let i = 0; i < atoms.length; i += 1) {
            if (selected.has(i)) {
              expect(shifted[i].x).toBeCloseTo(atoms[i].x + delta.x)
              expect(shifted[i].y).toBeCloseTo(atoms[i].y + delta.y)
              expect(shifted[i].z).toBeCloseTo(atoms[i].z + delta.z)
            } else {
              expect(shifted[i]).toEqual(atoms[i])
            }
          }
        },
      ),
      { numRuns: 120 },
    )
  })
})
