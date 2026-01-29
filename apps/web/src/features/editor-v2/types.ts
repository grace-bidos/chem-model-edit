import type { QeParameters, Structure } from '@/lib/types'

export type WorkspaceFile = {
  id: string
  name: string
  kind: 'in' | 'out'
  label: string
  pdbText?: string
  bcifUrl?: string
  structureId?: string
  structure?: Structure
  qeParams?: QeParameters | null
  parseSource?: string
  qeInput?: string
  initialOpenSections: {
    table: boolean
    parameter: boolean
  }
}

export type ToolMode = 'transfer' | 'supercell' | 'vibration'
