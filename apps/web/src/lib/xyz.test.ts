import { describe, expect, it } from 'vitest'

import { atomsToXyz, parseXyzBlock } from './xyz'

describe('xyz utilities', () => {
  it('ignores invalid and empty lines while parsing', () => {
    const parsed = parseXyzBlock('\nC 0 1 2\ninvalid\nH 3 4 5\n')

    expect(parsed).toHaveLength(2)
    expect(parsed[0]).toMatchObject({ symbol: 'C', x: 0, y: 1, z: 2 })
    expect(parsed[1]).toMatchObject({ symbol: 'H', x: 3, y: 4, z: 5 })
  })

  it('serializes atoms into xyz lines with fixed precision', () => {
    const xyz = atomsToXyz([
      { symbol: 'C', x: 1, y: 2.1234567, z: 3 },
      { symbol: 'H', x: 0, y: 0, z: 0.01 },
    ])

    expect(xyz).toContain('C 1.000000 2.123457 3.000000')
    expect(xyz).toContain('H 0.000000 0.000000 0.010000')
  })
})
