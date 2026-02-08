/* @vitest-environment jsdom */
import { cleanup, fireEvent, render, screen } from '@testing-library/react'
import { afterEach, describe, expect, it, vi } from 'vitest'

import { AtomTable } from './AtomTable'

afterEach(() => {
  cleanup()
})

describe('AtomTable', () => {
  it('renders headers and formatted coordinates', () => {
    render(
      <AtomTable
        rows={[
          { index: 0, symbol: 'Si', x: 1.23456, y: 2.3, z: 3.4 },
          { index: 1, symbol: 'O', x: 4.1, y: 5.2, z: 6.3 },
        ]}
        digits={2}
      />,
    )

    expect(screen.getByText('#')).toBeTruthy()
    expect(screen.getByText('El')).toBeTruthy()
    expect(screen.getByText('X')).toBeTruthy()
    expect(screen.getByText('1.23')).toBeTruthy()
    expect(screen.getByText('2.30')).toBeTruthy()
  })

  it('renders empty text when no rows are provided', () => {
    render(<AtomTable rows={[]} emptyText="No data" />)

    expect(screen.getByText('No data')).toBeTruthy()
  })

  it('applies selected and fixed row styles and click behavior', () => {
    const onRowClick = vi.fn<(index: number) => void>()

    render(
      <AtomTable
        rows={[
          { index: 0, symbol: 'Si', x: 0, y: 0, z: 0 },
          { index: 1, symbol: 'O', x: 1, y: 1, z: 1 },
        ]}
        selectedIndices={new Set([1])}
        fixedIndices={new Set([0])}
        onRowClick={onRowClick}
      />,
    )

    const fixedCell = screen.getByText('Si')
    const selectedCell = screen.getByText('O')
    const fixedRow = fixedCell.closest('tr')
    const selectedRow = selectedCell.closest('tr')

    expect(fixedRow?.className).toContain('bg-slate-50')
    expect(selectedRow?.className).toContain('bg-sky-100')

    fireEvent.click(fixedCell)
    fireEvent.click(selectedCell)

    expect(onRowClick).toHaveBeenCalledTimes(1)
    expect(onRowClick).toHaveBeenCalledWith(1)
  })
})
