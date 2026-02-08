import type { QeParameters, Structure } from '@/lib/types'

/** エディタワークスペースで管理するファイル状態。 */
export type WorkspaceFile = {
  id: string
  name: string
  kind: 'in' | 'out'
  label: string
  cifUrl?: string
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

/** ツールパネルの表示モード。 */
export type ToolMode = 'transfer' | 'supercell' | 'vibration'
