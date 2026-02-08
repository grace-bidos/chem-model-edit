/* @vitest-environment jsdom */
import { render, screen } from '@testing-library/react'
import { describe, expect, it } from 'vitest'

import { FilePanel } from './FilePanel'

import type { WorkspaceFile } from '../types'
import type { QeParameters, Structure } from '@/lib/types'

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
})
