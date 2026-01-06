import { createFileRoute } from '@tanstack/react-router'
import { useMemo, useState } from 'react'
import { Boxes, Grid3X3, Layers, Wand2 } from 'lucide-react'
import MolstarViewer from '../components/molstar/MolstarViewer'
import {
  generateTiledSupercell,
  latticeParamsToVectors,
  latticeVectorsToParams,
} from '../lib/api'
import { atomsToPdb } from '../lib/pdb'
import type { Atom, Lattice, LatticeParams, SupercellMeta } from '../lib/types'
import { parseXyzBlock } from '../lib/xyz'

export const Route = createFileRoute('/supercell')({
  component: SupercellPage,
})

function SupercellPage() {
  const sampleXyzA = ['O 0 0 0', 'H 0.757 0.586 0', 'H -0.757 0.586 0'].join(
    '\n',
  )
  const sampleXyzB = ['C 0 0 0', 'H 0.629 0.629 0.629'].join('\n')
  const [xyzA, setXyzA] = useState(sampleXyzA)
  const [xyzB, setXyzB] = useState(sampleXyzB)
  const [patternText, setPatternText] = useState('AB\nBA')
  const [patternGrid, setPatternGrid] = useState<string[][]>([
    ['A', 'B'],
    ['B', 'A'],
  ])
  const [patternError, setPatternError] = useState<string | null>(null)
  const [overlapCheck, setOverlapCheck] = useState(false)
  const [overlapTolerance, setOverlapTolerance] = useState('0.2')
  const [lattice, setLattice] = useState({
    a: { x: '5.0', y: '0.0', z: '0.0' },
    b: { x: '0.0', y: '5.0', z: '0.0' },
    c: { x: '0.0', y: '0.0', z: '5.0' },
  })
  const [atoms, setAtoms] = useState<Atom[]>([])
  const [meta, setMeta] = useState<SupercellMeta>({
    na: 0,
    nb: 0,
    layers: 0,
    overlapCount: 0,
  })
  const [error, setError] = useState<string | null>(null)
  const [converterError, setConverterError] = useState<string | null>(null)
  const [isGenerating, setIsGenerating] = useState(false)
  const [unit, setUnit] = useState<'angstrom' | 'bohr' | 'alat'>('angstrom')
  const [paramsDraft, setParamsDraft] = useState({
    a: '5.0',
    b: '5.0',
    c: '5.0',
    alpha: '90.0',
    beta: '90.0',
    gamma: '90.0',
  })
  const [vectorsDraft, setVectorsDraft] = useState({
    a: { x: '5.0', y: '0.0', z: '0.0' },
    b: { x: '0.0', y: '5.0', z: '0.0' },
    c: { x: '0.0', y: '0.0', z: '5.0' },
  })

  const BOHR_TO_ANGSTROM = 0.529177210903
  const convertUnit = (value: number, from: string, to: string) => {
    if (from === to) {
      return value
    }
    if (from === 'bohr' && to === 'angstrom') {
      return value * BOHR_TO_ANGSTROM
    }
    if (from === 'angstrom' && to === 'bohr') {
      return value / BOHR_TO_ANGSTROM
    }
    return value
  }

  const pdb = useMemo(() => atomsToPdb(atoms), [atoms])

  const gridToText = (grid: string[][]) =>
    grid.map((row) => row.join('')).join('\n')

  const parsePatternText = (text: string) => {
    const rows = text
      .split('\n')
      .map((line) => line.trim())
      .filter((line) => line.length > 0)
      .map((line) => line.replace(/\s+/g, '').toUpperCase())
    if (rows.length === 0) {
      return { grid: [] as string[][], error: 'タイルパターンが空です。' }
    }
    const grid = rows.map((row) => row.split(''))
    if (grid.some((row) => row.length === 0)) {
      return { grid: [] as string[][], error: 'タイルパターンが空です。' }
    }
    if (grid.some((row) => row.some((cell) => cell !== 'A' && cell !== 'B'))) {
      return { grid: [] as string[][], error: 'タイルは A/B のみで指定してください。' }
    }
    const width = grid[0].length
    if (grid.some((row) => row.length !== width)) {
      return { grid: [] as string[][], error: '行の長さが揃っていません。' }
    }
    return { grid, error: null }
  }

  const handlePatternTextChange = (value: string) => {
    setPatternText(value)
    const parsed = parsePatternText(value)
    if (parsed.error) {
      setPatternError(parsed.error)
      return
    }
    setPatternError(null)
    setPatternGrid(parsed.grid)
  }

  const updatePatternGrid = (next: string[][]) => {
    setPatternGrid(next)
    setPatternText(gridToText(next))
    setPatternError(null)
  }

  const handleToggleCell = (rowIndex: number, colIndex: number) => {
    updatePatternGrid(
      patternGrid.map((row, r) =>
        row.map((cell, c) => {
          if (r !== rowIndex || c !== colIndex) {
            return cell
          }
          return cell === 'A' ? 'B' : 'A'
        }),
      ),
    )
  }

  const handleAddRow = () => {
    if (patternGrid.length === 0) {
      updatePatternGrid([['A']])
      return
    }
    const width = patternGrid[0].length
    updatePatternGrid([...patternGrid, Array.from({ length: width }, () => 'A')])
  }

  const handleAddColumn = () => {
    if (patternGrid.length === 0) {
      updatePatternGrid([['A']])
      return
    }
    updatePatternGrid(patternGrid.map((row) => [...row, 'A']))
  }

  const handleUnitChange = (nextUnit: 'angstrom' | 'bohr' | 'alat') => {
    if (unit === nextUnit) {
      return
    }
    const scaleValue = (value: string) => {
      const parsed = Number(value)
      if (Number.isNaN(parsed)) {
        return value
      }
      return String(convertUnit(parsed, unit, nextUnit))
    }
    setParamsDraft((prev) => ({
      a: scaleValue(prev.a),
      b: scaleValue(prev.b),
      c: scaleValue(prev.c),
      alpha: prev.alpha,
      beta: prev.beta,
      gamma: prev.gamma,
    }))
    setVectorsDraft((prev) => ({
      a: { x: scaleValue(prev.a.x), y: scaleValue(prev.a.y), z: scaleValue(prev.a.z) },
      b: { x: scaleValue(prev.b.x), y: scaleValue(prev.b.y), z: scaleValue(prev.b.z) },
      c: { x: scaleValue(prev.c.x), y: scaleValue(prev.c.y), z: scaleValue(prev.c.z) },
    }))
    setUnit(nextUnit)
  }

  const handleGenerate = async () => {
    setIsGenerating(true)
    setError(null)
    try {
      const structureA = { atoms: parseXyzBlock(xyzA) }
      const structureB = { atoms: parseXyzBlock(xyzB) }
      const { a, b, c } = parseLatticeDraft(lattice)
      if (structureA.atoms.length === 0 || structureB.atoms.length === 0) {
        throw new Error('A/B両方の原子が必要です。')
      }
      if (patternError) {
        throw new Error(patternError)
      }
      if (patternGrid.length === 0) {
        throw new Error('タイルパターンが空です。')
      }
      const tolerance = overlapCheck ? Number(overlapTolerance) : undefined
      if (overlapCheck && (tolerance === undefined || Number.isNaN(tolerance))) {
        throw new Error('重複チェックの許容誤差が無効です。')
      }
      const result = await generateTiledSupercell({
        structureA,
        structureB,
        pattern: patternGrid,
        lattice: { a, b, c },
        checkOverlap: overlapCheck,
        overlapTolerance: tolerance,
      })
      setAtoms(result.structure.atoms)
      setMeta(result.meta)
    } catch (err) {
      setError(err instanceof Error ? err.message : '生成に失敗しました。')
    } finally {
      setIsGenerating(false)
    }
  }

  const parseLatticeDraft = (draft: typeof lattice): Lattice => {
    const parseVector = (vec: { x: string; y: string; z: string }) => ({
      x: Number(vec.x),
      y: Number(vec.y),
      z: Number(vec.z),
    })
    const a = parseVector(draft.a)
    const b = parseVector(draft.b)
    const c = parseVector(draft.c)
    if (
      [a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z].some((v) =>
        Number.isNaN(v),
      )
    ) {
      throw new Error('格子ベクトルの数値が無効です。')
    }
    return { a, b, c }
  }

  const parseParamsDraft = (): LatticeParams => {
    const values = {
      a: Number(paramsDraft.a),
      b: Number(paramsDraft.b),
      c: Number(paramsDraft.c),
      alpha: Number(paramsDraft.alpha),
      beta: Number(paramsDraft.beta),
      gamma: Number(paramsDraft.gamma),
    }
    if (Object.values(values).some((v) => Number.isNaN(v))) {
      throw new Error('格子定数の数値が無効です。')
    }
    return values
  }

  const handleVectorsToParams = async () => {
    setConverterError(null)
    try {
      const latticeVectors = parseLatticeDraft(vectorsDraft)
      const result = await latticeVectorsToParams({ lattice: latticeVectors, unit })
      setParamsDraft({
        a: result.params.a.toFixed(6),
        b: result.params.b.toFixed(6),
        c: result.params.c.toFixed(6),
        alpha: result.params.alpha.toFixed(4),
        beta: result.params.beta.toFixed(4),
        gamma: result.params.gamma.toFixed(4),
      })
    } catch (err) {
      setConverterError(err instanceof Error ? err.message : '変換に失敗しました。')
    }
  }

  const handleParamsToVectors = async () => {
    setConverterError(null)
    try {
      const params = parseParamsDraft()
      const result = await latticeParamsToVectors({ params, unit })
      setVectorsDraft({
        a: {
          x: result.lattice.a.x.toFixed(6),
          y: result.lattice.a.y.toFixed(6),
          z: result.lattice.a.z.toFixed(6),
        },
        b: {
          x: result.lattice.b.x.toFixed(6),
          y: result.lattice.b.y.toFixed(6),
          z: result.lattice.b.z.toFixed(6),
        },
        c: {
          x: result.lattice.c.x.toFixed(6),
          y: result.lattice.c.y.toFixed(6),
          z: result.lattice.c.z.toFixed(6),
        },
      })
    } catch (err) {
      setConverterError(err instanceof Error ? err.message : '変換に失敗しました。')
    }
  }

  return (
    <div className="min-h-screen bg-slate-950 text-white">
      <div className="mx-auto max-w-[1200px] px-6 py-10">
        <div className="flex flex-wrap items-center justify-between gap-4">
          <div>
            <p className="text-xs uppercase tracking-[0.3em] text-white/50">
              Supercell Builder
            </p>
            <h1 className="mt-2 text-3xl font-semibold">Stack slabs with intent</h1>
            <p className="mt-2 text-white/60">
              Define lattice vectors and tile patterns to generate a composite
              supercell preview.
            </p>
          </div>
          <button
            className="flex items-center gap-2 rounded-full bg-amber-300 px-4 py-2 text-sm font-semibold text-slate-900 shadow-lg shadow-amber-400/30"
            onClick={handleGenerate}
            disabled={isGenerating}
          >
            <Wand2 className="h-4 w-4" />
            {isGenerating ? 'Generating...' : 'Generate'}
          </button>
        </div>

        <div className="mt-8 grid gap-6 lg:grid-cols-[1.2fr_0.8fr]">
          <div className="rounded-2xl border border-white/10 bg-white/5 p-6">
            <div className="flex items-center gap-3 text-white/70">
              <Grid3X3 className="h-5 w-5" />
              <span className="text-sm uppercase tracking-[0.3em]">Lattice</span>
            </div>
            <div className="mt-6 grid grid-cols-3 gap-3 text-sm">
              {(['a', 'b', 'c'] as const).map((axis) => (
                <div
                  key={axis}
                  className="rounded-xl border border-white/10 bg-slate-900/60 px-3 py-4"
                >
                  <p className="text-xs text-white/40">{axis}-vector</p>
                  <div className="mt-2 space-y-1 text-xs text-white/60">
                    {(['x', 'y', 'z'] as const).map((component) => (
                      <input
                        key={`${axis}-${component}`}
                        value={lattice[axis][component]}
                        onChange={(event) =>
                          setLattice((prev) => ({
                            ...prev,
                            [axis]: {
                              ...prev[axis],
                              [component]: event.target.value,
                            },
                          }))
                        }
                        className="w-full rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                      />
                    ))}
                  </div>
                </div>
              ))}
            </div>
            <div className="mt-6 rounded-xl border border-white/10 bg-slate-900/60 p-4">
              <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                Tile Pattern
              </p>
              <textarea
                value={patternText}
                onChange={(event) => handlePatternTextChange(event.target.value)}
                className="mt-3 h-20 w-full resize-none rounded-md border border-white/10 bg-slate-950/70 p-2 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                aria-label="Tile pattern"
              />
              {patternError ? (
                <p className="mt-2 text-xs text-rose-300">{patternError}</p>
              ) : null}
              <div className="mt-4 flex flex-wrap items-center justify-between gap-2 text-xs text-white/60">
                <span>
                  {patternGrid.length} × {patternGrid[0]?.length ?? 0}
                </span>
                <div className="flex items-center gap-2">
                  <button
                    className="rounded-full border border-white/10 px-3 py-1 text-xs text-white/70 transition hover:text-white"
                    onClick={handleAddRow}
                    type="button"
                  >
                    + Row
                  </button>
                  <button
                    className="rounded-full border border-white/10 px-3 py-1 text-xs text-white/70 transition hover:text-white"
                    onClick={handleAddColumn}
                    type="button"
                  >
                    + Col
                  </button>
                </div>
              </div>
              <div className="mt-3 grid gap-2">
                {patternGrid.map((row, rowIndex) => (
                  <div key={`row-${rowIndex}`} className="flex gap-2">
                    {row.map((cell, colIndex) => (
                      <button
                        key={`cell-${rowIndex}-${colIndex}`}
                        className={`h-8 w-8 rounded-lg border text-xs font-semibold transition ${
                          cell === 'A'
                            ? 'border-sky-300/60 bg-sky-400/20 text-sky-100'
                            : 'border-rose-300/60 bg-rose-400/20 text-rose-100'
                        }`}
                        onClick={() => handleToggleCell(rowIndex, colIndex)}
                        type="button"
                      >
                        {cell}
                      </button>
                    ))}
                  </div>
                ))}
              </div>
              <div className="mt-4 flex items-center justify-between text-xs text-white/70">
                <label className="flex items-center gap-2">
                  <input
                    type="checkbox"
                    className="h-4 w-4 accent-amber-300"
                    checked={overlapCheck}
                    onChange={(event) => setOverlapCheck(event.target.checked)}
                  />
                  Overlap Check
                </label>
                <input
                  value={overlapTolerance}
                  onChange={(event) => setOverlapTolerance(event.target.value)}
                  className="w-20 rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  placeholder="tol"
                  disabled={!overlapCheck}
                />
              </div>
            </div>
            <div className="mt-6 grid gap-4 lg:grid-cols-2">
              <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Structure A (XYZ)
                </p>
                <textarea
                  value={xyzA}
                  onChange={(event) => setXyzA(event.target.value)}
                  className="mt-3 h-28 w-full resize-none rounded-md border border-white/10 bg-slate-950/70 p-2 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  aria-label="Structure A XYZ"
                />
              </div>
              <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
                <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                  Structure B (XYZ)
                </p>
                <textarea
                  value={xyzB}
                  onChange={(event) => setXyzB(event.target.value)}
                  className="mt-3 h-28 w-full resize-none rounded-md border border-white/10 bg-slate-950/70 p-2 font-mono text-xs text-white/80 focus:border-amber-300 focus:outline-none"
                  aria-label="Structure B XYZ"
                />
              </div>
            </div>
            {error ? <p className="mt-4 text-xs text-rose-300">{error}</p> : null}
          </div>

          <div className="flex flex-col gap-6">
            <div className="rounded-2xl border border-white/10 bg-white/5 p-6">
              <div className="flex items-center gap-3 text-white/70">
                <Boxes className="h-5 w-5" />
                <span className="text-sm uppercase tracking-[0.3em]">Output</span>
              </div>
              <div className="mt-6 space-y-3 text-sm text-white/70">
                <div className="flex items-center justify-between">
                  <span>Tile Grid</span>
                  <span className="font-semibold text-white">
                    {meta.na} × {meta.nb} × 1
                  </span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Atoms</span>
                  <span className="font-semibold text-white">
                    {atoms.length}
                  </span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Tiles</span>
                  <span className="font-semibold text-white">{meta.layers}</span>
                </div>
                <div className="flex items-center justify-between">
                  <span>Overlaps</span>
                  <span className="font-semibold text-white">
                    {meta.overlapCount ?? 0}
                  </span>
                </div>
              </div>
            </div>
            <div className="flex flex-1 flex-col rounded-2xl border border-white/10 bg-white/5 p-6">
              <div className="flex items-center justify-between text-white/60">
                <span className="text-xs uppercase tracking-[0.3em]">Viewer</span>
                <Layers className="h-4 w-4" />
              </div>
              <div className="mt-4 min-h-[240px] flex-1">
                <MolstarViewer
                  structures={
                    atoms.length
                      ? [{ id: 'supercell', pdbText: pdb, opacity: 1, visible: true }]
                      : []
                  }
                />
              </div>
            </div>
          </div>
        </div>

        <div className="mt-8 rounded-2xl border border-white/10 bg-white/5 p-6">
          <div className="flex items-center justify-between">
            <div>
              <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                Lattice Converter
              </p>
              <p className="mt-2 text-lg font-semibold">a,b,c & angles ⇄ vectors</p>
            </div>
            <select
              value={unit}
              onChange={(event) =>
                handleUnitChange(event.target.value as 'angstrom' | 'bohr' | 'alat')
              }
              className="rounded-full border border-white/10 bg-slate-950/70 px-3 py-1 text-xs text-white/80"
            >
              <option value="angstrom">angstrom</option>
              <option value="bohr">bohr</option>
              <option value="alat">alat</option>
            </select>
          </div>
          <div className="mt-6 grid gap-6 lg:grid-cols-2">
            <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
              <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                Lattice Params
              </p>
              <div className="mt-3 grid grid-cols-3 gap-2 text-xs text-white/70">
                {(['a', 'b', 'c'] as const).map((key) => (
                  <input
                    key={`param-${key}`}
                    value={paramsDraft[key]}
                    onChange={(event) =>
                      setParamsDraft((prev) => ({ ...prev, [key]: event.target.value }))
                    }
                    className="rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                  />
                ))}
                {(['alpha', 'beta', 'gamma'] as const).map((key) => (
                  <input
                    key={`param-${key}`}
                    value={paramsDraft[key]}
                    onChange={(event) =>
                      setParamsDraft((prev) => ({ ...prev, [key]: event.target.value }))
                    }
                    className="rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                  />
                ))}
              </div>
              <button
                className="mt-4 w-full rounded-xl border border-white/10 px-3 py-2 text-xs text-white/70 transition hover:text-white"
                onClick={handleParamsToVectors}
              >
                Params → Vectors
              </button>
            </div>
            <div className="rounded-xl border border-white/10 bg-slate-900/60 p-4">
              <p className="text-xs uppercase tracking-[0.3em] text-white/50">
                Cell Vectors
              </p>
              <div className="mt-3 space-y-2 text-xs text-white/70">
                {(['a', 'b', 'c'] as const).map((axis) => (
                  <div key={`vec-${axis}`} className="grid grid-cols-3 gap-2">
                    {(['x', 'y', 'z'] as const).map((component) => (
                      <input
                        key={`${axis}-${component}`}
                        value={vectorsDraft[axis][component]}
                        onChange={(event) =>
                          setVectorsDraft((prev) => ({
                            ...prev,
                            [axis]: {
                              ...prev[axis],
                              [component]: event.target.value,
                            },
                          }))
                        }
                        className="rounded-md border border-white/10 bg-slate-950/70 px-2 py-1 text-white/80 focus:border-amber-300 focus:outline-none"
                      />
                    ))}
                  </div>
                ))}
              </div>
              <button
                className="mt-4 w-full rounded-xl border border-amber-200/40 px-3 py-2 text-xs text-amber-100/90 transition hover:text-amber-200"
                onClick={handleVectorsToParams}
              >
                Vectors → Params
              </button>
            </div>
          </div>
          {converterError ? (
            <p className="mt-4 text-xs text-rose-300">{converterError}</p>
          ) : null}
        </div>
      </div>
    </div>
  )
}
