import type { QeParameters, Structure } from '@/lib/types'

import type { WorkspaceFile } from '../types'

export const mockStructureA: Structure = {
  atoms: [
    { symbol: 'Si', x: 0, y: 0, z: 0 },
    { symbol: 'O', x: 1.6, y: 0, z: 0 },
    { symbol: 'O', x: -1.6, y: 0, z: 0 },
  ],
  lattice: {
    a: { x: 5.43, y: 0, z: 0 },
    b: { x: 0, y: 5.43, z: 0 },
    c: { x: 0, y: 0, z: 5.43 },
  },
}

export const mockStructureB: Structure = {
  atoms: [
    { symbol: 'C', x: 0, y: 0, z: 0 },
    { symbol: 'H', x: 0.9, y: 0.9, z: 0 },
    { symbol: 'H', x: -0.9, y: 0.9, z: 0 },
    { symbol: 'H', x: 0.9, y: -0.9, z: 0 },
    { symbol: 'H', x: -0.9, y: -0.9, z: 0 },
  ],
  lattice: {
    a: { x: 8, y: 0, z: 0 },
    b: { x: 0, y: 8, z: 0 },
    c: { x: 0, y: 0, z: 8 },
  },
}

export const mockQeParams: QeParameters = {
  control: { calculation: 'relax', pseudo_dir: './pseudo' },
  system: { ecutwfc: 40, ecutrho: 320 },
  electrons: { conv_thr: 1e-8 },
  kpoints: { grid: [3, 3, 1], shift: [0, 0, 0] },
}

export const emptyWorkspaceFile: WorkspaceFile = {
  id: 'ws-empty',
  name: 'empty.in',
  kind: 'in',
  label: 'Empty Input',
  initialOpenSections: {
    table: true,
    parameter: true,
  },
}

export const populatedWorkspaceFile: WorkspaceFile = {
  id: 'ws-a',
  name: 'si-o.in',
  kind: 'in',
  label: 'SiO slab',
  structureId: 'structure-sio',
  structure: mockStructureA,
  qeParams: mockQeParams,
  initialOpenSections: {
    table: true,
    parameter: true,
  },
}

export const mockWorkspaceFiles: Array<WorkspaceFile> = [
  populatedWorkspaceFile,
  {
    id: 'ws-b',
    name: 'ch4.in',
    kind: 'in',
    label: 'CH4 molecule',
    structureId: 'structure-ch4',
    structure: mockStructureB,
    qeParams: null,
    initialOpenSections: {
      table: false,
      parameter: false,
    },
  },
]
