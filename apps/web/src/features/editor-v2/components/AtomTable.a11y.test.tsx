/* @vitest-environment jsdom */
import { render } from '@testing-library/react'
import { axe } from 'jest-axe'
import { describe, expect, it } from 'vitest'

import { AtomTable } from './AtomTable'

describe('AtomTable a11y', () => {
  it('has no basic accessibility violations', async () => {
    const { container } = render(
      <AtomTable
        rows={[
          { index: 0, symbol: 'Si', x: 0.1234, y: 1.2345, z: 2.3456 },
          { index: 1, symbol: 'O', x: 3.4567, y: 4.5678, z: 5.6789 },
        ]}
        digits={3}
      />, 
    )
    const results = await axe(container)
    expect(results.violations).toEqual([])
  })
})
