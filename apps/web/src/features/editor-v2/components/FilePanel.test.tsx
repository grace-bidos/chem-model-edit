/* @vitest-environment jsdom */
import { cleanup, fireEvent, render, screen, waitFor } from '@testing-library/react'
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

import { FilePanel } from './FilePanel'

import type { WorkspaceFile } from '../types'
import type { QeParameters, Structure } from '@/lib/types'

const { mockedGetStructure } = vi.hoisted(() => ({
  mockedGetStructure: vi.fn(),
}))

vi.mock('@/lib/api', () => ({
  getStructure: mockedGetStructure,
}))

vi.mock('@/components/molstar/MolstarViewer', () => ({
  default: ({
    onError,
    onLoad,
  }: {
    onError?: (message: string) => void
    onLoad?: () => void
  }) => (
    <div>
      <button type="button" onClick={() => onLoad?.()}>
        trigger-viewer-load
      </button>
      <button type="button" onClick={() => onError?.('viewer boom')}>
        trigger-viewer-error
      </button>
    </div>
  ),
}))

const baseStructure = {
  atoms: [{ symbol: 'Si', x: 0, y: 0, z: 0 }],
} as unknown as Structure

const buildFile = (overrides: Partial<WorkspaceFile> = {}): WorkspaceFile => ({
  id: 'file-1',
  name: 'sample.in',
  kind: 'in',
  label: 'Sample',
  structure: baseStructure,
  qeParams: null,
  initialOpenSections: {
    table: true,
    parameter: true,
  },
  ...overrides,
})

describe('FilePanel parameter table', () => {
  beforeEach(() => {
    vi.clearAllMocks()
  })

  afterEach(() => {
    cleanup()
  })

  it('renders grouped parameter tables and formatted values', () => {
    const params = {
      control: {
        calculation: 'scf',
      },
      system: {
        ecutwfc: 40,
        flags: [1, 2, 3],
        metadata: { mode: 'test' },
      },
      electrons: {},
      ions: {},
      cell: {},
      pseudopotentials: {},
      kpoints: { mode: 'gamma' },
    } as unknown as QeParameters

    render(<FilePanel data={buildFile({ qeParams: params })} fileId="file-1" />)

    expect(screen.getByText('Control')).toBeTruthy()
    expect(screen.getByText('System')).toBeTruthy()
    expect(screen.getByText('K-points')).toBeTruthy()
    expect(screen.getByText('scf')).toBeTruthy()
    expect(screen.getByText('1, 2, 3')).toBeTruthy()
    expect(screen.getByText('{"mode":"test"}')).toBeTruthy()
  })

  it('shows empty message when parameters are not available', () => {
    render(<FilePanel data={buildFile()} fileId="file-1" />)

    expect(screen.getByText('No parameters.')).toBeTruthy()
  })

  it('loads structure/params when structureId is present and data is missing', async () => {
    const loadedStructure = {
      atoms: [{ symbol: 'O', x: 1.2, y: 2.3, z: 3.4 }],
    } as unknown as Structure
    const loadedParams = {
      control: { calculation: 'relax' },
      system: {},
      electrons: {},
      ions: {},
      cell: {},
      pseudopotentials: {},
      kpoints: {},
    } as unknown as QeParameters
    mockedGetStructure.mockResolvedValueOnce({
      structure: loadedStructure,
      params: loadedParams,
    })
    const onStructureLoaded = vi.fn<(fileId: string, structure: Structure) => void>()
    const onParamsLoaded = vi.fn<
      (fileId: string, params: QeParameters | null) => void
    >()

    render(
      <FilePanel
        data={buildFile({
          structure: undefined,
          qeParams: null,
          structureId: 'structure-1',
        })}
        fileId="file-1"
        onStructureLoaded={onStructureLoaded}
        onParamsLoaded={onParamsLoaded}
      />,
    )

    await waitFor(() => expect(mockedGetStructure).toHaveBeenCalledWith('structure-1'))
    await waitFor(() => expect(screen.getByText('O')).toBeTruthy())
    expect(screen.getByText('Control')).toBeTruthy()
    expect(screen.getByText('relax')).toBeTruthy()
    expect(onStructureLoaded).toHaveBeenCalledWith('file-1', loadedStructure)
    expect(onParamsLoaded).toHaveBeenCalledWith('file-1', loadedParams)
  })

  it('shows load error states when fetching structure fails', async () => {
    mockedGetStructure.mockRejectedValueOnce(new Error('network down'))

    render(
      <FilePanel
        data={buildFile({
          structure: undefined,
          qeParams: null,
          structureId: 'structure-err',
        })}
        fileId="file-1"
      />,
    )

    await waitFor(() => expect(mockedGetStructure).toHaveBeenCalled())
    await waitFor(() =>
      expect(screen.getByText('Failed to load structure.')).toBeTruthy(),
    )
    expect(screen.getByText('Failed to load parameters.')).toBeTruthy()
  })

  it('does not fetch when structure and params are already present', () => {
    render(
      <FilePanel
        data={buildFile({
          structureId: 'already-hydrated',
          qeParams: {
            control: { calculation: 'scf' },
            system: {},
            electrons: {},
            ions: {},
            cell: {},
            pseudopotentials: {},
            kpoints: {},
          } as unknown as QeParameters,
        })}
        fileId="file-1"
      />,
    )

    expect(mockedGetStructure).not.toHaveBeenCalled()
  })

  it('renders header controls and triggers callbacks', () => {
    const onClose = vi.fn()
    const onMinimize = vi.fn()
    render(
      <FilePanel
        data={buildFile()}
        fileId="file-1"
        onClose={onClose}
        onMinimize={onMinimize}
      />,
    )

    fireEvent.click(screen.getByLabelText('Minimize file panel'))
    fireEvent.click(screen.getByLabelText('Close file panel'))
    expect(onMinimize).toHaveBeenCalledTimes(1)
    expect(onClose).toHaveBeenCalledTimes(1)
  })

  it('shows and clears viewer error overlay for cif viewer', async () => {
    render(
      <FilePanel
        data={buildFile({
          cifUrl: 'https://example.test/structure.cif',
        })}
        fileId="file-1"
      />,
    )

    fireEvent.click(screen.getByRole('button', { name: 'trigger-viewer-error' }))
    expect(screen.getByText('Viewer failed to load')).toBeTruthy()
    expect(screen.getByText('viewer boom')).toBeTruthy()

    fireEvent.click(screen.getByRole('button', { name: 'trigger-viewer-load' }))
    await waitFor(() =>
      expect(screen.queryByText('Viewer failed to load')).toBeNull(),
    )
  })
})
