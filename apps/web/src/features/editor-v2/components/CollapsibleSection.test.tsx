/* @vitest-environment jsdom */
import { fireEvent, render, screen } from '@testing-library/react'
import { describe, expect, it } from 'vitest'

import { CollapsibleSection } from './CollapsibleSection'

describe('CollapsibleSection', () => {
  it('respects defaultOpen and custom classes', () => {
    const { container } = render(
      <CollapsibleSection
        title="Parameters"
        defaultOpen
        className="outer-class"
        contentClassName="content-class"
      >
        <div>body content</div>
      </CollapsibleSection>,
    )

    expect(screen.getByText('body content')).toBeTruthy()
    expect(container.firstElementChild?.className).toContain('outer-class')
    expect(container.querySelector('.content-class')).toBeTruthy()
  })

  it('toggles open state on header click', () => {
    render(
      <CollapsibleSection title="Table">
        <div>rows</div>
      </CollapsibleSection>,
    )

    expect(screen.queryByText('rows')).toBeNull()

    fireEvent.click(screen.getByRole('button', { name: 'Table' }))
    expect(screen.getByText('rows')).toBeTruthy()

    fireEvent.click(screen.getByRole('button', { name: 'Table' }))
    expect(screen.queryByText('rows')).toBeNull()
  })
})
