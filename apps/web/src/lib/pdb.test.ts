import { describe, expect, it } from 'vitest'

import { atomsToPdb } from './pdb'

describe('atomsToPdb', () => {
  it('formats atoms into a PDB-like block', () => {
    const pdb = atomsToPdb([
      { symbol: 'H', x: 0, y: 0, z: 0 },
      { symbol: 'O', x: 0.5, y: 0.5, z: 0.5 },
    ])
    expect(pdb).toContain('HEADER')
    expect(pdb).toContain('ATOM')
    expect(pdb).toContain('H')
    expect(pdb).toContain('O')
    expect(pdb).toContain('END')
  })
})
