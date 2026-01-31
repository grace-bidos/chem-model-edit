import type {
  AuthMe,
  AuthSession,
  DeltaTransplantResponse,
  Lattice,
  LatticeParams,
  Structure,
  StructureCreateResponse,
  StructureGetResponse,
  StructureParseResponse,
  StructureExportResponse,
  SupercellBuildRequest,
  SupercellBuildResponse,
  SupercellMeta,
  SupercellResponse,
  ZPEJobRequest,
  ZPEJobResponse,
  ZPEJobStatus,
  ZPEParseResponse,
  ZPEQueueTargetList,
  ZPEResult,
} from './types'

export type ApiRequest = {
  path: string
  method?: 'GET' | 'POST' | 'PUT' | 'PATCH' | 'DELETE'
  body?: unknown
  token?: string | null
  responseType?: 'json' | 'text'
}

export type ApiRequester = <T>(options: ApiRequest) => Promise<T>

export type ApiClientOptions = {
  request: ApiRequester
  getToken?: () => string | undefined
}

export type ApiClient = {
  parseQeInput: (content: string) => Promise<Structure>
  createStructureFromQe: (content: string) => Promise<StructureCreateResponse>
  getStructure: (structureId: string) => Promise<StructureGetResponse>
  exportQeInput: (structure: Structure) => Promise<string>
  deltaTransplant: (params: { smallIn: string; smallOut: string; largeIn: string }) => Promise<string>
  generateSupercell: (params: {
    structureA: Structure
    structureB: Structure
    sequence: string
    lattice: Lattice
  }) => Promise<{ structure: Structure; meta: SupercellMeta }>
  generateTiledSupercell: (params: {
    structureA: Structure
    structureB: Structure
    pattern: Array<Array<string>>
    lattice: Lattice
    checkOverlap: boolean
    overlapTolerance?: number
  }) => Promise<{ structure: Structure; meta: SupercellMeta }>
  buildSupercell: (params: SupercellBuildRequest) => Promise<SupercellBuildResponse>
  latticeVectorsToParams: (params: { lattice: Lattice; unit: string }) => Promise<{ lattice: Lattice; params: LatticeParams; unit: string }>
  latticeParamsToVectors: (params: { params: LatticeParams; unit: string }) => Promise<{ lattice: Lattice; params: LatticeParams; unit: string }>
  parseZpeInput: (content: string, structureId?: string | null) => Promise<ZPEParseResponse>
  registerAccount: (params: { email: string; password: string }) => Promise<AuthSession>
  loginAccount: (params: { email: string; password: string }) => Promise<AuthSession>
  fetchAuthMe: () => Promise<AuthMe>
  logoutAccount: () => Promise<void>
  createEnrollToken: (params?: { ttlSeconds?: number; label?: string }) => Promise<{ token: string; expiresAt: string; ttlSeconds: number; label?: string | null }>
  fetchQueueTargets: (params?: { limit?: number; offset?: number }) => Promise<ZPEQueueTargetList>
  selectQueueTarget: (targetId: string) => Promise<{ activeTargetId: string }>
  createZpeJob: (request: ZPEJobRequest) => Promise<ZPEJobResponse>
  fetchZpeStatus: (jobId: string) => Promise<ZPEJobStatus>
  fetchZpeResult: (jobId: string) => Promise<ZPEResult>
  downloadZpeFile: (jobId: string, kind: 'summary' | 'freqs') => Promise<string>
}

export const createApiClient = (options: ApiClientOptions): ApiClient => {
  const { request, getToken } = options
  const withToken = (token?: string | null) => token ?? getToken?.()

  const requestJson = async <T>(params: ApiRequest): Promise<T> => {
    return request<T>({ ...params, responseType: 'json' })
  }

  const requestText = async (params: ApiRequest): Promise<string> => {
    return request<string>({ ...params, responseType: 'text' })
  }

  return {
    async parseQeInput(content) {
      const data = await requestJson<StructureParseResponse>({
        path: '/structures/parse',
        method: 'POST',
        body: { content },
      })
      return data.structure
    },

    async createStructureFromQe(content) {
      return requestJson<StructureCreateResponse>({
        path: '/structures',
        method: 'POST',
        body: { content },
      })
    },

    async getStructure(structureId) {
      const safeId = encodeURIComponent(structureId)
      return requestJson<StructureGetResponse>({
        path: `/structures/${safeId}`,
      })
    },

    async exportQeInput(structure) {
      const data = await requestJson<StructureExportResponse>({
        path: '/structures/export',
        method: 'POST',
        body: { structure },
      })
      return data.content
    },

    async deltaTransplant(params) {
      const data = await requestJson<DeltaTransplantResponse>({
        path: '/transforms/delta-transplant',
        method: 'POST',
        body: {
          smallIn: params.smallIn,
          smallOut: params.smallOut,
          largeIn: params.largeIn,
        },
      })
      return data.content
    },

    async generateSupercell(params) {
      return requestJson<SupercellResponse>({
        path: '/supercells',
        method: 'POST',
        body: params,
      })
    },

    async generateTiledSupercell(params) {
      return requestJson<SupercellResponse>({
        path: '/supercells/tiled',
        method: 'POST',
        body: params,
      })
    },

    async buildSupercell(params) {
      return requestJson<SupercellBuildResponse>({
        path: '/supercells/builds',
        method: 'POST',
        body: params,
      })
    },

    async latticeVectorsToParams(params) {
      return requestJson<{ lattice: Lattice; params: LatticeParams; unit: string }>({
        path: '/lattices/convert',
        method: 'POST',
        body: { from: 'vectors', lattice: params.lattice, unit: params.unit },
      })
    },

    async latticeParamsToVectors(params) {
      return requestJson<{ lattice: Lattice; params: LatticeParams; unit: string }>({
        path: '/lattices/convert',
        method: 'POST',
        body: { from: 'params', params: params.params, unit: params.unit },
      })
    },

    async parseZpeInput(content, structureId) {
      const payload: { content: string; structureId?: string } = { content }
      if (structureId) {
        payload.structureId = structureId
      }
      return requestJson<ZPEParseResponse>({
        path: '/zpe/parse',
        method: 'POST',
        body: payload,
      })
    },

    async registerAccount(params) {
      return requestJson<AuthSession>({
        path: '/auth/register',
        method: 'POST',
        body: params,
      })
    },

    async loginAccount(params) {
      return requestJson<AuthSession>({
        path: '/auth/login',
        method: 'POST',
        body: params,
      })
    },

    async fetchAuthMe() {
      return requestJson<AuthMe>({
        path: '/auth/me',
        token: withToken(),
      })
    },

    async logoutAccount() {
      await requestJson({
        path: '/auth/logout',
        method: 'POST',
        token: withToken(),
      })
    },

    async createEnrollToken(params) {
      return requestJson({
        path: '/zpe/compute/enroll-tokens',
        method: 'POST',
        token: withToken(),
        body: {
          ttlSeconds: params?.ttlSeconds,
          label: params?.label,
        },
      })
    },

    async fetchQueueTargets(params) {
      const search = new URLSearchParams()
      if (params?.limit !== undefined) {
        search.set('limit', String(params.limit))
      }
      if (params?.offset !== undefined) {
        search.set('offset', String(params.offset))
      }
      const query = search.toString()
      return requestJson<ZPEQueueTargetList>({
        path: `/zpe/targets${query ? `?${query}` : ''}`,
        token: withToken(),
      })
    },

    async selectQueueTarget(targetId) {
      const safeId = encodeURIComponent(targetId)
      return requestJson<{ activeTargetId: string }>({
        path: `/zpe/targets/${safeId}/active`,
        method: 'PUT',
        token: withToken(),
      })
    },

    async createZpeJob(requestPayload) {
      return requestJson<ZPEJobResponse>({
        path: '/zpe/jobs',
        method: 'POST',
        token: withToken(),
        body: requestPayload,
      })
    },

    async fetchZpeStatus(jobId) {
      const safeId = encodeURIComponent(jobId)
      return requestJson<ZPEJobStatus>({
        path: `/zpe/jobs/${safeId}`,
        token: withToken(),
      })
    },

    async fetchZpeResult(jobId) {
      const safeId = encodeURIComponent(jobId)
      const data = await requestJson<{ result: ZPEResult }>({
        path: `/zpe/jobs/${safeId}/result`,
        token: withToken(),
      })
      return data.result
    },

    async downloadZpeFile(jobId, kind) {
      const safeId = encodeURIComponent(jobId)
      return requestText({
        path: `/zpe/jobs/${safeId}/files?kind=${kind}`,
        token: withToken(),
      })
    },
  }
}
