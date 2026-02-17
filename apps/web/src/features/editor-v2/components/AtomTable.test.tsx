/* @vitest-environment jsdom */
import { cleanup, fireEvent, render, screen } from '@testing-library/react'
import { afterEach, describe, expect, it, vi } from 'vitest'

import { AtomTable } from './AtomTable'

afterEach(() => {
  cleanup()
})

describe('AtomTable', () => {
  const rows = [
    { index: 0, symbol: 'Si', x: 1.23456, y: 2.3, z: 3.4 },
    { index: 1, symbol: 'O', x: 4.1, y: 5.2, z: 6.3 },
  ]

  it('renders headers and formatted coordinates', () => {
    render(<AtomTable rows={rows} digits={2} />)

    expect(screen.getByText('#')).toBeTruthy()
    expect(screen.getByText('El')).toBeTruthy()
    expect(screen.getByText('X')).toBeTruthy()
    expect(screen.getByText('1.23')).toBeTruthy()
    expect(screen.getByText('2.30')).toBeTruthy()
  })

  it('renders default empty text when no rows are provided', () => {
    render(<AtomTable rows={[]} />)
    expect(screen.getByText('No atoms.')).toBeTruthy()
  })

  it('respects showIndex and stickyHeader toggles', () => {
    const { rerender, container } = render(
      <AtomTable rows={rows} showIndex={false} stickyHeader={false} />,
    )

    expect(screen.queryByText('#')).toBeNull()
    const header = container.querySelector('thead')
    expect(header?.className).not.toContain('sticky')

    rerender(<AtomTable rows={rows} showIndex stickyHeader />)

    expect(screen.getByText('#')).toBeTruthy()
    expect(header?.className).toContain('sticky')
  })

  it('applies selected and fixed row styles and click behavior', () => {
    const onRowClick = vi.fn<(index: number) => void>()

    render(
      <AtomTable
        rows={rows}
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
    expect(selectedRow?.className).toContain('cursor-pointer')

    fireEvent.click(fixedCell)
    fireEvent.click(selectedCell)

    expect(onRowClick).toHaveBeenCalledTimes(1)
    expect(onRowClick).toHaveBeenCalledWith(1)
  })

  it('uses row.fixed as fallback and disables clicks for fixed rows', () => {
    const onRowClick = vi.fn<(index: number) => void>()
    render(
      <AtomTable
        rows={[{ index: 0, symbol: 'Fe', x: 0, y: 0, z: 0, fixed: true }]}
        onRowClick={onRowClick}
      />,
    )

    const symbolCell = screen.getByText('Fe')
    const row = symbolCell.closest('tr')
    const coordCells = screen.getAllByText('0.000')

    expect(row?.className).toContain('bg-slate-50')
    expect(row?.className).not.toContain('cursor-pointer')
    expect(coordCells[0].className).toContain('text-slate-500')

    fireEvent.click(symbolCell)
    expect(onRowClick).not.toHaveBeenCalled()
  })

  it('does not allow selection when selectionEnabled is false', () => {
    const onRowClick = vi.fn<(index: number) => void>()

    render(<AtomTable rows={rows} onRowClick={onRowClick} selectionEnabled={false} />)

    const selectedCell = screen.getByText('O')
    const selectedRow = selectedCell.closest('tr')
    expect(selectedRow?.className).not.toContain('cursor-pointer')

    fireEvent.click(selectedCell)
    expect(onRowClick).not.toHaveBeenCalled()
  })
})
