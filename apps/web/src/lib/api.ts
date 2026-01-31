import { getAuthToken } from './auth'

import type {
  AuthMe,
  AuthSession,
  Lattice,
  LatticeParams,
  QeParameters,
  Structure,
  SupercellBuildRequest,
  SupercellBuildResponse,
  SupercellMeta,
  ZPEJobRequest,
  ZPEJobResponse,
  ZPEJobStatus,
  ZPEParseResponse,
  ZPEQueueTargetList,
  ZPEResult,
} from './types'

import { requestApi } from '@/server/api'

const API_BASE = import.meta.env.VITE_API_BASE ?? 'http://localhost:8000'

type RequestOptions = Parameters<typeof requestApi>[0]['data']
type RequestInput = Omit<RequestOptions, 'responseType'>

const requestJson = async <T>(params: RequestInput): Promise<T> => {
  return (await requestApi({
    data: {
      ...params,
      responseType: 'json',
    },
  })) as T
}

const requestText = async (params: RequestInput): Promise<string> => {
  return (await requestApi({
    data: {
      ...params,
      responseType: 'text',
    },
  })) as string
}

const getToken = () => getAuthToken() ?? undefined

export async function parseQeInput(content: string): Promise<Structure> {
  const data = await requestJson<{ structure: Structure }>({
    path: '/parse',
    method: 'POST',
    body: { content },
  })
  return data.structure
}

export async function createStructureFromQe(content: string): Promise<{
  structure_id: string
  structure: Structure
  source: string
  params?: QeParameters | null
  raw_input?: string | null
}> {
  return requestJson<{
    structure_id: string
    structure: Structure
    source: string
    params?: QeParameters | null
    raw_input?: string | null
  }>({
    path: '/structures',
    method: 'POST',
    body: { content },
  })
}

export async function getStructure(structureId: string): Promise<{
  structure: Structure
  params?: QeParameters | null
  raw_input?: string | null
  source?: string | null
}> {
  const safeId = encodeURIComponent(structureId)
  return requestJson<{
    structure: Structure
    params?: QeParameters | null
    raw_input?: string | null
    source?: string | null
  }>({
    path: `/structures/${safeId}`,
  })
}

export function structureViewUrl(
  structureId: string,
  params?: {
    format?: 'cif'
  },
): string {
  const format = params?.format ?? 'cif'
  const query = new URLSearchParams({
    format,
  })
  const safeId = encodeURIComponent(structureId)
  return `${API_BASE}/structures/${safeId}/view?${query.toString()}`
}

export async function exportQeInput(structure: Structure): Promise<string> {
  const data = await requestJson<{ content: string }>({
    path: '/export',
    method: 'POST',
    body: { structure },
  })
  return data.content
}

export async function deltaTransplant(params: {
  smallIn: string
  smallOut: string
  largeIn: string
}): Promise<string> {
  const data = await requestJson<{ content: string }>({
    path: '/transplant/delta',
    method: 'POST',
    body: {
      small_in: params.smallIn,
      small_out: params.smallOut,
      large_in: params.largeIn,
    },
  })
  return data.content
}

export async function generateSupercell(params: {
  structureA: Structure
  structureB: Structure
  sequence: string
  lattice: Lattice
}): Promise<{ structure: Structure; meta: SupercellMeta }> {
  return requestJson<{ structure: Structure; meta: SupercellMeta }>({
    path: '/supercell',
    method: 'POST',
    body: params,
  })
}

export async function generateTiledSupercell(params: {
  structureA: Structure
  structureB: Structure
  pattern: Array<Array<string>>
  lattice: Lattice
  checkOverlap: boolean
  overlapTolerance?: number
}): Promise<{ structure: Structure; meta: SupercellMeta }> {
  return requestJson<{ structure: Structure; meta: SupercellMeta }>({
    path: '/supercell/tiled',
    method: 'POST',
    body: params,
  })
}

export async function buildSupercell(
  params: SupercellBuildRequest,
): Promise<SupercellBuildResponse> {
  return requestJson<SupercellBuildResponse>({
    path: '/supercell/build',
    method: 'POST',
    body: params,
  })
}

export async function latticeVectorsToParams(params: {
  lattice: Lattice
  unit: string
}): Promise<{ lattice: Lattice; params: LatticeParams; unit: string }> {
  return requestJson<{
    lattice: Lattice
    params: LatticeParams
    unit: string
  }>({
    path: '/lattice/vectors-to-params',
    method: 'POST',
    body: params,
  })
}

export async function latticeParamsToVectors(params: {
  params: LatticeParams
  unit: string
}): Promise<{ lattice: Lattice; params: LatticeParams; unit: string }> {
  return requestJson<{
    lattice: Lattice
    params: LatticeParams
    unit: string
  }>({
    path: '/lattice/params-to-vectors',
    method: 'POST',
    body: params,
  })
}

export async function parseZpeInput(
  content: string,
  structureId?: string | null,
): Promise<ZPEParseResponse> {
  const payload: { content: string; structure_id?: string } = { content }
  if (structureId) {
    payload.structure_id = structureId
  }
  return requestJson<ZPEParseResponse>({
    path: '/calc/zpe/parse',
    method: 'POST',
    body: payload,
  })
}

export async function registerAccount(params: {
  email: string
  password: string
}): Promise<AuthSession> {
  return requestJson<AuthSession>({
    path: '/auth/register',
    method: 'POST',
    body: params,
  })
}

export async function loginAccount(params: {
  email: string
  password: string
}): Promise<AuthSession> {
  return requestJson<AuthSession>({
    path: '/auth/login',
    method: 'POST',
    body: params,
  })
}

export async function fetchAuthMe(): Promise<AuthMe> {
  return requestJson<AuthMe>({
    path: '/auth/me',
    token: getToken(),
  })
}

export async function logoutAccount(): Promise<void> {
  await requestJson<{ ok: boolean }>({
    path: '/auth/logout',
    method: 'POST',
    token: getToken(),
  })
}

export async function createEnrollToken(params?: {
  ttlSeconds?: number
  label?: string
}): Promise<{
  token: string
  expires_at: string
  ttl_seconds: number
  label?: string | null
}> {
  return requestJson({
    path: '/calc/zpe/compute/enroll-tokens',
    method: 'POST',
    token: getToken(),
    body: {
      ttl_seconds: params?.ttlSeconds,
      label: params?.label,
    },
  })
}

export async function fetchQueueTargets(): Promise<ZPEQueueTargetList> {
  return requestJson<ZPEQueueTargetList>({
    path: '/calc/zpe/compute/targets',
    token: getToken(),
  })
}

export async function selectQueueTarget(
  targetId: string,
): Promise<{ active_target_id: string }> {
  return requestJson({
    path: '/calc/zpe/compute/targets/select',
    method: 'POST',
    token: getToken(),
    body: { target_id: targetId },
  })
}

export async function createZpeJob(
  request: ZPEJobRequest,
): Promise<ZPEJobResponse> {
  return requestJson<ZPEJobResponse>({
    path: '/calc/zpe/jobs',
    method: 'POST',
    token: getToken(),
    body: request,
  })
}

export async function fetchZpeStatus(jobId: string): Promise<ZPEJobStatus> {
  const safeId = encodeURIComponent(jobId)
  return requestJson<ZPEJobStatus>({
    path: `/calc/zpe/jobs/${safeId}`,
    token: getToken(),
  })
}

export async function fetchZpeResult(jobId: string): Promise<ZPEResult> {
  const safeId = encodeURIComponent(jobId)
  const data = await requestJson<{ result: ZPEResult }>({
    path: `/calc/zpe/jobs/${safeId}/result`,
    token: getToken(),
  })
  return data.result
}

export async function downloadZpeFile(
  jobId: string,
  kind: 'summary' | 'freqs',
): Promise<string> {
  const safeId = encodeURIComponent(jobId)
  return requestText({
    path: `/calc/zpe/jobs/${safeId}/files?kind=${kind}`,
    token: getToken(),
  })
}
