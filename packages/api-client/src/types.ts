import type { components, paths } from './generated/schema'

type Schemas = components['schemas']

export type OpenApiPaths = paths
export type OpenApiSchemas = Schemas

export type Atom = Schemas['Atom']
export type Vector3 = Schemas['Vector3']
export type Lattice = Schemas['Lattice']
export type Structure = Schemas['Structure-Output']
export type StructureInput = Schemas['Structure-Input']
export type QeParameters = Schemas['QeParameters']
export type Pagination = Schemas['Pagination']

export type StructureParseResponse = Schemas['StructureParseResponse']
export type StructureCreateResponse = Schemas['StructureCreateResponse']
export type StructureGetResponse = Schemas['StructureGetResponse']
export type StructureExportResponse = Schemas['StructureExportResponse']
export type StructureParseRequest = Schemas['StructureParseRequest']
export type StructureCreateRequest = Schemas['StructureCreateRequest']
export type StructureExportRequest = Schemas['StructureExportRequest']

export type DeltaTransplantRequest = Schemas['DeltaTransplantRequest']
export type DeltaTransplantResponse = Schemas['DeltaTransplantResponse']

export type SupercellGridAxis = Schemas['SupercellGridAxis']
export type SupercellGrid = Schemas['SupercellGrid']
export type SupercellBuildOptions = Schemas['SupercellBuildOptions']
export type SupercellBuildOutput = Schemas['SupercellBuildOutput']
export type SupercellBuildRequest = Schemas['SupercellBuildRequest']
export type SupercellBuildMeta = Schemas['SupercellBuildMeta']
export type SupercellBuildResponse = Schemas['SupercellBuildResponse']

export type ZPEParseResponse = Schemas['ZPEParseResponse']
export type ZPEJobStatus = Schemas['ZPEJobStatus']
export type ZPEQueueTarget = Schemas['QueueTarget']
export type ZPEQueueTargetList = Schemas['QueueTargetListResponse']
export type ZPEQueueTargetSelectResponse = Schemas['QueueTargetSelectResponse']

export type ExecutionEvent = Schemas['ExecutionEvent']
export type RuntimeEventAck = Schemas['RuntimeEventAck']
export type RuntimeJobStatus = Schemas['RuntimeJobStatusResponse']
export type SubmitJobAccepted = Schemas['SubmitJobAccepted']
export type SubmitJobCommand = Schemas['SubmitJobCommand']

export type RuntimeJobDetail = {
  [key: string]: unknown
}

// Backward-compatibility type for web UI rendering.
// Runtime-only contract no longer exposes this schema directly.
export type ZPEResult = {
  calc_type?: string
  freqs_cm: number[]
  zpe_ev: number
  s_vib_jmol_k: number
  mobile_indices: number[]
  fixed_indices: number[]
  kpoints: [number, number, number]
  delta: number
  low_cut_cm: number
  temperature: number
  use_environ: boolean
  qe_input: string
  pseudo_dir: string
  calc_start_time: string
  calc_end_time: string
  elapsed_seconds: number
  cache_checked: number
  cache_deleted: number
  ecutwfc?: number | null
  ecutrho?: number | null
}
