import type {
  AuthLoginRequest,
  AuthMe,
  AuthRegisterRequest,
  AuthSession,
  DeltaTransplantRequest,
  DeltaTransplantResponse,
  EnrollTokenRequest,
  EnrollTokenResponse,
  Structure,
  StructureCreateResponse,
  StructureExportRequest,
  StructureGetResponse,
  StructureParseRequest,
  SupercellBuildRequest,
  SupercellBuildResponse,
  ZPEJobRequest,
  ZPEJobResponse,
  ZPEJobStatus,
  ZPEParseResponse,
  ZPEQueueTargetList,
  ZPEQueueTargetSelectResponse,
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
  getStructure: (structure_id: string) => Promise<StructureGetResponse>
  exportQeInput: (structure: Structure) => Promise<string>
  deltaTransplant: (params: DeltaTransplantRequest) => Promise<string>
  buildSupercell: (params: SupercellBuildRequest) => Promise<SupercellBuildResponse>
  parseZpeInput: (content: string, structure_id?: string | null) => Promise<ZPEParseResponse>
  registerAccount: (params: AuthRegisterRequest) => Promise<AuthSession>
  loginAccount: (params: AuthLoginRequest) => Promise<AuthSession>
  fetchAuthMe: () => Promise<AuthMe>
  logoutAccount: () => Promise<void>
  createEnrollToken: (params?: EnrollTokenRequest) => Promise<EnrollTokenResponse>
  fetchQueueTargets: (params?: { limit?: number; offset?: number }) => Promise<ZPEQueueTargetList>
  selectQueueTarget: (target_id: string) => Promise<ZPEQueueTargetSelectResponse>
  createZpeJob: (request: ZPEJobRequest) => Promise<ZPEJobResponse>
  fetchZpeStatus: (job_id: string) => Promise<ZPEJobStatus>
  fetchZpeResult: (job_id: string) => Promise<ZPEResult>
  downloadZpeFile: (job_id: string, kind: 'summary' | 'freqs') => Promise<string>
  structureViewPath: (
    structure_id: string,
    params?: { format?: 'cif' },
  ) => string
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
      const body: StructureParseRequest = { content }
      const data = await requestJson<{ structure: Structure }>({
        path: '/structures/parse',
        method: 'POST',
        body,
      })
      return data.structure
    },

    async createStructureFromQe(content) {
      const body: StructureParseRequest = { content }
      return requestJson<StructureCreateResponse>({
        path: '/structures',
        method: 'POST',
        body,
      })
    },

    async getStructure(structure_id) {
      const safeId = encodeURIComponent(structure_id)
      return requestJson<StructureGetResponse>({
        path: `/structures/${safeId}`,
      })
    },

    async exportQeInput(structure) {
      const body: StructureExportRequest = { structure }
      const data = await requestJson<{ content: string }>({
        path: '/structures/export',
        method: 'POST',
        body,
      })
      return data.content
    },

    async deltaTransplant(params) {
      const data = await requestJson<DeltaTransplantResponse>({
        path: '/transforms/delta-transplant',
        method: 'POST',
        body: params,
      })
      return data.content
    },

    async buildSupercell(params) {
      return requestJson<SupercellBuildResponse>({
        path: '/supercells/builds',
        method: 'POST',
        body: params,
      })
    },

    async parseZpeInput(content, structure_id) {
      const payload: { content: string; structure_id?: string } = { content }
      if (structure_id) {
        payload.structure_id = structure_id
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
      return requestJson<EnrollTokenResponse>({
        path: '/zpe/compute/enroll-tokens',
        method: 'POST',
        token: withToken(),
        body: {
          ttl_seconds: params?.ttl_seconds,
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

    async selectQueueTarget(target_id) {
      const safeId = encodeURIComponent(target_id)
      return requestJson<ZPEQueueTargetSelectResponse>({
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

    async fetchZpeStatus(job_id) {
      const safeId = encodeURIComponent(job_id)
      return requestJson<ZPEJobStatus>({
        path: `/zpe/jobs/${safeId}`,
        token: withToken(),
      })
    },

    async fetchZpeResult(job_id) {
      const safeId = encodeURIComponent(job_id)
      const data = await requestJson<{ result: ZPEResult }>({
        path: `/zpe/jobs/${safeId}/result`,
        token: withToken(),
      })
      return data.result
    },

    async downloadZpeFile(job_id, kind) {
      const safeId = encodeURIComponent(job_id)
      return requestText({
        path: `/zpe/jobs/${safeId}/files?kind=${kind}`,
        token: withToken(),
      })
    },

    structureViewPath(structure_id, params) {
      const format = params?.format ?? 'cif'
      const safeId = encodeURIComponent(structure_id)
      const query = new URLSearchParams({ format })
      return `/structures/${safeId}/view?${query.toString()}`
    },
  }
}
