import fc from 'fast-check'
import { describe, expect, it } from 'vitest'

import { atomsToXyz, parseXyzBlock } from './xyz'

import type { Atom } from './types'

const finiteFloat = fc.float({ min: -1000, max: 1000, noNaN: true })
const symbol = fc.stringMatching(/^[A-Z][a-z]{0,2}$/)

const atomArb = fc.record({
  symbol,
  x: finiteFloat,
  y: finiteFloat,
  z: finiteFloat,
}) satisfies fc.Arbitrary<Atom>

describe('atomsToXyz/parseXyzBlock property tests', () => {
  it('preserves symbol and six-decimal coordinates through round-trip', () => {
    fc.assert(
      fc.property(fc.array(atomArb, { minLength: 1, maxLength: 20 }), (atoms) => {
        const serialized = atomsToXyz(atoms)
        const parsed = parseXyzBlock(serialized)

        expect(parsed).toHaveLength(atoms.length)
        for (let i = 0; i < atoms.length; i += 1) {
          expect(parsed[i].symbol).toBe(atoms[i].symbol)
          expect(parsed[i].x).toBe(Number(atoms[i].x.toFixed(6)))
          expect(parsed[i].y).toBe(Number(atoms[i].y.toFixed(6)))
          expect(parsed[i].z).toBe(Number(atoms[i].z.toFixed(6)))
        }
      }),
      { numRuns: 150 },
    )
  })
})
