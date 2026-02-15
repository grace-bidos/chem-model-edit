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
  TableContainer,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'

/** {@link AtomTable} が表示する1行分の原子データ。 */
export type AtomTableRow = {
  index: number
  symbol: string
  x: number
  y: number
  z: number
  fixed?: boolean
}

/** {@link AtomTable} が受け取る props。 */
export type AtomTableProps = {
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

const columnHelper = createColumnHelper<AtomTableRow>()
const coordinateColumns = new Set(['x', 'y', 'z'])

const formatNumber = (value: number, digits: number) => value.toFixed(digits)

const buildColumns = (
  digits: number,
  showIndex: boolean,
): Array<ColumnDef<AtomTableRow>> => {
  const columns: Array<ColumnDef<AtomTableRow>> = []

  if (showIndex) {
    columns.push(
      columnHelper.accessor('index', {
        header: '#',
        cell: ({ getValue }) => getValue() + 1,
      }) as ColumnDef<AtomTableRow>,
    )
  }

  columns.push(
    columnHelper.accessor('symbol', {
      header: 'El',
    }) as ColumnDef<AtomTableRow>,
  )

  columns.push(
    columnHelper.accessor('x', {
      header: 'X',
      cell: ({ getValue }) => formatNumber(getValue(), digits),
    }) as ColumnDef<AtomTableRow>,
    columnHelper.accessor('y', {
      header: 'Y',
      cell: ({ getValue }) => formatNumber(getValue(), digits),
    }) as ColumnDef<AtomTableRow>,
    columnHelper.accessor('z', {
      header: 'Z',
      cell: ({ getValue }) => formatNumber(getValue(), digits),
    }) as ColumnDef<AtomTableRow>,
  )

  return columns
}

type RowVisualState = {
  clickable: boolean
  isFixed: boolean
  isSelected: boolean
}

const getRowVisualState = (
  row: AtomTableRow,
  canSelect: boolean,
  selectedIndices?: Set<number>,
  fixedIndices?: Set<number>,
): RowVisualState => {
  const isFixed = fixedIndices?.has(row.index) ?? row.fixed ?? false
  const isSelected = selectedIndices?.has(row.index) ?? false
  return {
    isFixed,
    isSelected,
    clickable: canSelect && !isFixed,
  }
}

/** TanStack の行モデルと shadcn Table で原子座標テーブルを描画する。 */
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
  const columns = useMemo(
    () => buildColumns(digits, showIndex),
    [digits, showIndex],
  )

  const table = useReactTable({
    data: rows,
    columns,
    getCoreRowModel: getCoreRowModel(),
  })

  const colSpan = table.getVisibleLeafColumns().length

  return (
    <TableContainer
      className={cn(
        'h-full rounded border border-slate-100 bg-slate-50 p-2',
        containerClassName,
      )}
    >
      <Table className="text-xs">
        <TableHeader
          className={cn(stickyHeader && 'sticky top-0 z-10 bg-slate-50')}
        >
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
              const visualState = getRowVisualState(
                row.original,
                canSelect,
                selectedIndices,
                fixedIndices,
              )
              const rowClasses = cn(
                visualState.isFixed && 'bg-slate-50',
                visualState.isSelected && 'bg-sky-100',
                visualState.clickable && 'cursor-pointer hover:bg-sky-200/40',
              )
              const fixedCoordClass = visualState.isFixed
                ? 'bg-slate-50 text-slate-500'
                : 'text-slate-600'

              return (
                <TableRow
                  key={row.id}
                  className={rowClasses}
                  onClick={
                    visualState.clickable && onRowClick
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
    </TableContainer>
  )
}
