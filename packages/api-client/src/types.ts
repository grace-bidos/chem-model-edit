export type Atom = {
  symbol: string
  x: number
  y: number
  z: number
}

export type Vector3 = {
  x: number
  y: number
  z: number
}

export type Lattice = {
  a: Vector3
  b: Vector3
  c: Vector3
}

export type LatticeParams = {
  a: number
  b: number
  c: number
  alpha: number
  beta: number
  gamma: number
}

export type Structure = {
  atoms: Atom[]
  lattice?: Lattice | null
}

export type QeParameters = {
  control: Record<string, unknown>
  system: Record<string, unknown>
  electrons: Record<string, unknown>
  ions: Record<string, unknown>
  cell: Record<string, unknown>
  pseudopotentials: Record<string, string>
  kpoints?: Record<string, unknown> | null
}

export type Pagination = {
  total: number
  limit: number
  offset: number
}

export type AuthUser = {
  id: string
  email: string
  createdAt: string
}

export type AuthSession = {
  token: string
  expiresAt: string
  user: AuthUser
}

export type AuthMe = {
  user: AuthUser
  expiresAt: string
}

export type AuthLogoutResponse = {
  ok: true
}

export type StructureParseResponse = {
  structure: Structure
}

export type StructureCreateResponse = {
  id: string
  structure: Structure
  source: string
  params?: QeParameters | null
  rawInput?: string | null
}

export type StructureGetResponse = {
  structure: Structure
  params?: QeParameters | null
  rawInput?: string | null
  source?: string | null
}

export type StructureExportResponse = {
  content: string
}

export type DeltaTransplantResponse = {
  content: string
}

export type SupercellMeta = {
  na: number
  nb: number
  layers: number
  overlapCount?: number
}

export type SupercellResponse = {
  structure: Structure
  meta: SupercellMeta
}

export type SupercellGridAxis =
  | { row: 'a'; col: 'b' }
  | { row: 'b'; col: 'a' }

export type SupercellGrid = {
  rows: number
  cols: number
  tiles: string[][]
  axis?: SupercellGridAxis
}

export type SupercellBuildOptions = {
  checkOverlap?: boolean
  overlapTolerance?: number
  validateLattice?: 'none' | 'warn' | 'error'
}

export type SupercellBuildOutput = {
  includeStructure?: boolean
}

export type SupercellBuildRequest = {
  baseStructureId: string
  grid: SupercellGrid
  options?: SupercellBuildOptions
  output?: SupercellBuildOutput
}

export type SupercellBuildMeta = {
  rows: number
  cols: number
  tileCount: number
  overlapCount?: number
  baseStructureId: string
  structureIdsUsed: string[]
}

export type SupercellBuildResponse = {
  id: string
  structure?: Structure
  meta: SupercellBuildMeta
}

export type ZPEParseResponse = {
  structure: Structure
  fixedIndices: number[]
  atomicSpecies: Record<string, string>
  kpoints?: [number, number, number] | null
}

export type ZPEJobRequest = {
  content: string
  mobileIndices: number[]
  useEnviron?: boolean
  inputDir?: string | null
  calcMode?: 'new' | 'continue'
  structureId?: string | null
}

export type ZPEJobResponse = {
  id: string
}

export type ZPEJobStatus = {
  status: string
  detail?: string | null
  updatedAt?: string | null
}

export type ZPEQueueTarget = {
  id: string
  queueName: string
  serverId: string
  registeredAt: string
  name?: string | null
}

export type ZPEQueueTargetList = {
  targets: ZPEQueueTarget[]
  activeTargetId?: string | null
  pagination: Pagination
}

export type ZPEResult = {
  freqsCm: number[]
  zpeEv: number
  sVibJmolK: number
  mobileIndices: number[]
  fixedIndices: number[]
  kpoints: [number, number, number]
  delta: number
  lowCutCm: number
  temperature: number
  useEnviron: boolean
  qeInput: string
  pseudoDir: string
  calcStartTime: string
  calcEndTime: string
  elapsedSeconds: number
  cacheChecked: number
  cacheDeleted: number
  ecutwfc?: number | null
  ecutrho?: number | null
}

export type ErrorResponse = {
  error: {
    code: string
    message: string
    details?: unknown
  }
}
