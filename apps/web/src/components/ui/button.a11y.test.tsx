/* @vitest-environment jsdom */
import { render } from '@testing-library/react'
import { axe } from 'jest-axe'
import { describe, expect, it } from 'vitest'

import { Button } from './button'

describe('Button a11y', () => {
  it('has no basic accessibility violations', async () => {
    const { container } = render(<Button>Run</Button>)
    const results = await axe(container)
    expect(results.violations).toEqual([])
  })
})
