import type { components, paths } from './generated/schema'

type Schemas = components['schemas']

export type OpenApiPaths = paths
export type OpenApiSchemas = Schemas

export type Atom = Schemas['Atom']
export type Vector3 = Schemas['Vector3']
export type Lattice = Schemas['Lattice']
export type LatticeParams = Schemas['LatticeParams']
export type Structure = Schemas['Structure-Output']
export type StructureInput = Schemas['Structure-Input']
export type QeParameters = Schemas['QeParameters']
export type Pagination = Schemas['Pagination']

export type AuthUser = Schemas['AuthUser']
export type AuthSession = Schemas['AuthSession']
export type AuthMe = Schemas['AuthMe']
export type AuthLogoutResponse = Schemas['AuthLogoutResponse']

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
export type ZPEJobRequest = Schemas['ZPEJobRequest']
export type ZPEJobResponse = Schemas['ZPEJobResponse']
export type ZPEJobStatus = Schemas['ZPEJobStatus']
export type ZPEQueueTarget = Schemas['QueueTarget']
export type ZPEQueueTargetList = Schemas['QueueTargetListResponse']
export type ZPEQueueTargetSelectResponse = Schemas['QueueTargetSelectResponse']
export type ZPEResult = Schemas['ZPEResult']
export type EnrollTokenResponse = Schemas['EnrollTokenResponse']
export type AuthLoginRequest = Schemas['AuthLoginRequest']
export type AuthRegisterRequest = Schemas['AuthRegisterRequest']
export type EnrollTokenRequest = Schemas['EnrollTokenRequest']
