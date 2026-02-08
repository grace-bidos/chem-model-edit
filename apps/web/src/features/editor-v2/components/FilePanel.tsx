import { useEffect, useMemo, useRef, useState } from 'react'
import { Cuboid, Minus, X } from 'lucide-react'

import { CollapsibleSection } from './CollapsibleSection'
import { AtomTable } from './AtomTable'

import type { WorkspaceFile } from '../types'
import type { QeParameters, Structure } from '@/lib/types'

import MolstarViewer from '@/components/molstar/MolstarViewer'
import {
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableHeader,
  TableRow,
} from '@/components/ui/table'
import { getStructure } from '@/lib/api'
import { cn } from '@/lib/utils'

interface FilePanelProps {
  data: WorkspaceFile
  fileId: string
  onStructureLoaded?: (fileId: string, structure: Structure) => void
  onParamsLoaded?: (fileId: string, params: QeParameters | null) => void
  onClose?: () => void
  onMinimize?: () => void
  showHeader?: boolean
  className?: string
}

type ParamGroup = {
  label: string
  entries: Array<[string, unknown]>
}

const formatParamValue = (value: unknown) => {
  if (value === null || value === undefined) {
    return ''
  }
  if (typeof value === 'number' || typeof value === 'boolean') {
    return String(value)
  }
  if (Array.isArray(value)) {
    return value.map((entry) => String(entry)).join(', ')
  }
  if (typeof value === 'object') {
    try {
      return JSON.stringify(value)
    } catch (_error) {
      return String(value)
    }
  }
  return String(value)
}

const buildParamGroups = (params: QeParameters | null): Array<ParamGroup> => {
  if (!params) {
    return []
  }
  const rawGroups: Array<{
    label: string
    values: Record<string, unknown> | null | undefined
  }> = [
    { label: 'Control', values: params.control },
    { label: 'System', values: params.system },
    { label: 'Electrons', values: params.electrons },
    { label: 'Ions', values: params.ions },
    { label: 'Cell', values: params.cell },
    { label: 'Pseudopotentials', values: params.pseudopotentials },
    { label: 'K-points', values: params.kpoints ?? {} },
  ]
  return rawGroups
    .map(({ label, values }) => ({
      label,
      entries: Object.entries(values ?? {}),
    }))
    .filter((group) => group.entries.length > 0)
}

type ParameterTableProps = {
  group: ParamGroup
}

/** Parameters セクション用のキー/値テーブルを描画する。 */
function ParameterTable({ group }: ParameterTableProps) {
  return (
    <TableContainer className="rounded border border-slate-100 bg-slate-50">
      <Table className="text-xs">
        <TableHeader className="sticky top-0 z-10 bg-slate-50">
          <TableRow>
            <TableHead className="w-2/5">Key</TableHead>
            <TableHead>Value</TableHead>
          </TableRow>
        </TableHeader>
        <TableBody>
          {group.entries.map(([key, value]) => (
            <TableRow key={`${group.label}-${key}`}>
              <TableCell className="font-medium text-slate-500">{key}</TableCell>
              <TableCell className="break-all font-mono text-[11px] text-slate-700">
                {formatParamValue(value) || '-'}
              </TableCell>
            </TableRow>
          ))}
        </TableBody>
      </Table>
    </TableContainer>
  )
}

/** ビューア・原子表・QEパラメータ表を含む単一ファイルパネルを表示する。 */
export function FilePanel({
  data,
  fileId,
  onStructureLoaded,
  onParamsLoaded,
  onClose,
  onMinimize,
  showHeader = true,
  className,
}: FilePanelProps) {
  const [viewerError, setViewerError] = useState<string | null>(null)
  const [structure, setStructure] = useState<Structure | null>(
    data.structure ?? null,
  )
  const [qeParams, setQeParams] = useState<QeParameters | null>(
    data.qeParams ?? null,
  )
  const [tableError, setTableError] = useState<string | null>(null)
  const [isTableLoading, setIsTableLoading] = useState(false)
  const [paramsError, setParamsError] = useState<string | null>(null)
  const [isParamsLoading, setIsParamsLoading] = useState(false)
  const onStructureLoadedRef =
    useRef<FilePanelProps['onStructureLoaded']>(onStructureLoaded)
  const onParamsLoadedRef =
    useRef<FilePanelProps['onParamsLoaded']>(onParamsLoaded)

  useEffect(() => {
    setViewerError(null)
  }, [data.cifUrl])

  useEffect(() => {
    onStructureLoadedRef.current = onStructureLoaded
  }, [onStructureLoaded])

  useEffect(() => {
    onParamsLoadedRef.current = onParamsLoaded
  }, [onParamsLoaded])

  useEffect(() => {
    setStructure(data.structure ?? null)
    setQeParams(data.qeParams ?? null)
  }, [data.structure, data.structureId, data.qeParams])

  useEffect(() => {
    const hasStructure = Boolean(data.structure)
    const hasParams = Boolean(data.qeParams)

    if (!data.structureId) {
      setIsTableLoading(false)
      setTableError(null)
      setIsParamsLoading(false)
      setParamsError(null)
      return
    }

    if (hasStructure) {
      setIsTableLoading(false)
      setTableError(null)
    }
    if (hasParams) {
      setIsParamsLoading(false)
      setParamsError(null)
    }

    const needsStructure = !hasStructure
    const needsParams = !hasParams
    if (!needsStructure && !needsParams) {
      return
    }

    let cancelled = false
    if (needsStructure) {
      setIsTableLoading(true)
      setTableError(null)
    }
    if (needsParams) {
      setIsParamsLoading(true)
      setParamsError(null)
    }

    getStructure(data.structureId)
      .then((result) => {
        if (cancelled) return
        if (needsStructure) {
          setStructure(result.structure)
          onStructureLoadedRef.current?.(fileId, result.structure)
        }
        if (needsParams) {
          setQeParams(result.params ?? null)
          onParamsLoadedRef.current?.(fileId, result.params ?? null)
        }
        if (needsStructure) {
          setIsTableLoading(false)
        }
        if (needsParams) {
          setIsParamsLoading(false)
        }
      })
      .catch((err) => {
        if (cancelled) return
        const message = err instanceof Error ? err.message : String(err)
        if (needsStructure) {
          setTableError(message)
          setIsTableLoading(false)
        }
        if (needsParams) {
          setParamsError(message)
          setIsParamsLoading(false)
        }
      })

    return () => {
      cancelled = true
    }
  }, [data.structure, data.structureId, data.qeParams, fileId])

  const atoms = structure?.atoms ?? []
  const atomRows = atoms.map((atom, index) => ({
    index,
    symbol: atom.symbol,
    x: atom.x,
    y: atom.y,
    z: atom.z,
  }))
  const tableEmptyText = isTableLoading
    ? 'Loading...'
    : tableError
      ? 'Failed to load structure.'
      : 'No atoms.'
  const paramsEmptyText = isParamsLoading
    ? 'Loading...'
    : paramsError
      ? 'Failed to load parameters.'
      : 'No parameters.'

  const paramGroups = useMemo(() => buildParamGroups(qeParams), [qeParams])

  return (
    <div
      className={cn(
        'flex h-full w-80 flex-shrink-0 flex-col gap-4 border-r border-border bg-background p-4 animate-in fade-in slide-in-from-left-4 duration-300',
        className,
      )}
    >
      {showHeader ? (
        <div className="flex items-center justify-between">
          <h2
            className="flex-1 truncate text-lg font-semibold tracking-tight"
            title={data.name}
          >
            {data.name}
          </h2>
          {onClose || onMinimize ? (
            <div className="flex items-center gap-1">
              {onMinimize ? (
                <button
                  type="button"
                  onClick={onMinimize}
                  className="rounded p-1 text-slate-400 transition-colors hover:bg-slate-100 hover:text-slate-600"
                  title="Minimize"
                  aria-label="Minimize file panel"
                >
                  <Minus className="h-4 w-4" aria-hidden="true" />
                </button>
              ) : null}
              {onClose ? (
                <button
                  type="button"
                  onClick={onClose}
                  className="rounded p-1 text-slate-400 transition-colors hover:bg-red-50 hover:text-red-500"
                  title="Close"
                  aria-label="Close file panel"
                >
                  <X className="h-4 w-4" aria-hidden="true" />
                </button>
              ) : null}
            </div>
          ) : null}
        </div>
      ) : null}

      <div className="flex min-h-0 flex-1 flex-col gap-4">
        <div className="group relative flex h-1/2 min-h-[220px] w-full flex-col overflow-hidden rounded-lg border border-border bg-card">
          <div className="absolute inset-0 bg-grid-slate-200/50 [mask-image:linear-gradient(0deg,white,rgba(255,255,255,0.7))]" />
          <div className="relative z-10 flex w-full flex-1 flex-col">
            {data.cifUrl ? (
              <MolstarViewer
                cifUrl={data.cifUrl}
                onError={setViewerError}
                onLoad={() => setViewerError(null)}
              />
            ) : (
              <div className="flex h-full w-full flex-col items-center justify-center px-4 text-center text-muted-foreground">
                <Cuboid className="mb-2 h-12 w-12 text-slate-300 transition-colors duration-300 group-hover:text-blue-400" />
                <span className="font-medium text-slate-600">{data.label}</span>
                <span className="text-xs text-slate-400">
                  3D Structure View
                </span>
              </div>
            )}
          </div>
          {viewerError ? (
            <div className="absolute inset-0 z-20 flex items-center justify-center bg-white/80 px-4 text-center">
              <div className="rounded-md border border-red-200 bg-red-50 px-3 py-2 text-sm text-red-700 shadow-sm">
                <p className="font-semibold">Viewer failed to load</p>
                <p className="mt-1 text-xs text-red-600">{viewerError}</p>
              </div>
            </div>
          ) : null}
        </div>

        <div className="flex min-h-0 flex-1 flex-col gap-2 overflow-hidden pr-1">
          <CollapsibleSection
            title="Table"
            defaultOpen={data.initialOpenSections.table}
            className="flex min-h-0 flex-1 flex-col"
            contentClassName="flex min-h-0 flex-1 flex-col"
          >
            <div className="flex min-h-0 flex-1 flex-col gap-2">
              <p className="text-xs text-muted-foreground">
                Atomic Positions (Angstrom)
              </p>
              <div className="min-h-0 flex-1">
                <AtomTable
                  rows={atomRows}
                  digits={4}
                  emptyText={tableEmptyText}
                  showIndex
                  containerClassName="h-full"
                />
              </div>
            </div>
          </CollapsibleSection>

          <div className="shrink-0">
            <CollapsibleSection
              title="Parameters"
              defaultOpen={data.initialOpenSections.parameter}
            >
              <div className="space-y-2">
                <p className="mb-1 text-xs text-muted-foreground">
                  QE Option Params
                </p>
                {paramGroups.length === 0 ? (
                  <div className="rounded border border-slate-100 bg-slate-50 px-3 py-2 text-xs text-slate-500">
                    {paramsEmptyText}
                  </div>
                ) : (
                  <div className="space-y-3">
                    {paramGroups.map((group) => (
                      <div key={group.label} className="space-y-2">
                        <p className="text-[11px] font-semibold uppercase tracking-wide text-slate-500">
                          {group.label}
                        </p>
                        <ParameterTable group={group} />
                      </div>
                    ))}
                  </div>
                )}
              </div>
            </CollapsibleSection>
          </div>
        </div>
      </div>
    </div>
  )
}
