import { useMemo } from 'react'
import {
  createColumnHelper,
  flexRender,
  getCoreRowModel,
  useReactTable,
} from '@tanstack/react-table'
import type { ColumnDef } from '@tanstack/react-table'

import { cn } from '@/lib/utils'
import {
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'

type AtomTableRow = {
  index: number
  symbol: string
  x: number
  y: number
  z: number
  fixed?: boolean
}

type AtomTableProps = {
  rows: Array<AtomTableRow>
  selectedIndices?: Set<number>
  fixedIndices?: Set<number>
  onRowClick?: (index: number) => void
  selectionEnabled?: boolean
  digits?: number
  emptyText?: string
  containerClassName?: string
  stickyHeader?: boolean
  showIndex?: boolean
}

const formatNumber = (value: number, digits: number) => value.toFixed(digits)
const columnHelper = createColumnHelper<AtomTableRow>()
const coordinateColumns = new Set(['x', 'y', 'z'])

export function AtomTable({
  rows,
  selectedIndices,
  fixedIndices,
  onRowClick,
  selectionEnabled = true,
  digits = 3,
  emptyText = 'No atoms.',
  containerClassName,
  stickyHeader = false,
  showIndex = true,
}: AtomTableProps) {
  const canSelect = Boolean(onRowClick) && selectionEnabled
  const columns = useMemo<Array<ColumnDef<AtomTableRow, any>>>(() => {
    const cols: Array<ColumnDef<AtomTableRow, any>> = []
    if (showIndex) {
      cols.push(
        columnHelper.accessor('index', {
          header: '#',
          cell: ({ getValue }) => getValue<number>() + 1,
        }),
      )
    }
    cols.push(
      columnHelper.accessor('symbol', {
        header: 'El',
      }),
    )
    cols.push(
      columnHelper.accessor('x', {
        header: 'X',
        cell: ({ getValue }) => formatNumber(getValue<number>(), digits),
      }),
      columnHelper.accessor('y', {
        header: 'Y',
        cell: ({ getValue }) => formatNumber(getValue<number>(), digits),
      }),
      columnHelper.accessor('z', {
        header: 'Z',
        cell: ({ getValue }) => formatNumber(getValue<number>(), digits),
      }),
    )
    return cols
  }, [digits, showIndex])

  const table = useReactTable({
    data: rows,
    columns,
    getCoreRowModel: getCoreRowModel(),
  })

  const colSpan = table.getVisibleLeafColumns().length
  return (
    <Table
      className="text-xs"
      containerClassName={cn('h-full', containerClassName)}
    >
      <TableHeader className={cn(stickyHeader && 'sticky top-0 bg-slate-50')}>
        {table.getHeaderGroups().map((headerGroup) => (
          <TableRow key={headerGroup.id}>
            {headerGroup.headers.map((header) => (
              <TableHead key={header.id}>
                {header.isPlaceholder
                  ? null
                  : flexRender(
                      header.column.columnDef.header,
                      header.getContext(),
                    )}
              </TableHead>
            ))}
          </TableRow>
        ))}
      </TableHeader>
      <TableBody>
        {table.getRowModel().rows.length ? (
          table.getRowModel().rows.map((row) => {
            const isFixed =
              fixedIndices?.has(row.original.index) ?? row.original.fixed ?? false
            const isSelected =
              selectedIndices?.has(row.original.index) ?? false
            const clickable = canSelect && !isFixed
            const rowClasses = cn(
              isFixed && 'bg-slate-50',
              isSelected && 'bg-sky-100',
              clickable && 'cursor-pointer hover:bg-sky-200/40',
            )
            const fixedCoordClass = isFixed
              ? 'bg-slate-50 text-slate-500'
              : 'text-slate-600'
            return (
              <TableRow
                key={row.id}
                className={rowClasses}
                onClick={
                  clickable && onRowClick
                    ? () => onRowClick(row.original.index)
                    : undefined
                }
              >
                {row.getVisibleCells().map((cell) => {
                  const columnId = cell.column.id
                  const cellClass = cn(
                    columnId === 'index' &&
                      'text-xs font-medium text-slate-500',
                    columnId === 'symbol' && 'font-semibold text-slate-700',
                    coordinateColumns.has(columnId) &&
                      cn('font-mono', fixedCoordClass),
                  )
                  return (
                    <TableCell key={cell.id} className={cellClass}>
                      {flexRender(
                        cell.column.columnDef.cell,
                        cell.getContext(),
                      )}
                    </TableCell>
                  )
                })}
              </TableRow>
            )
          })
        ) : (
          <TableRow>
            <TableCell
              colSpan={colSpan}
              className="py-10 text-center text-sm text-slate-500"
            >
              {emptyText}
            </TableCell>
          </TableRow>
        )}
      </TableBody>
    </Table>
  )
}
