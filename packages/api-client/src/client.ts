import type {
  DeltaTransplantRequest,
  DeltaTransplantResponse,
  ExecutionEvent,
  RuntimeEventAck,
  RuntimeJobDetail,
  RuntimeJobStatus,
  Structure,
  StructureCreateResponse,
  StructureExportRequest,
  StructureGetResponse,
  StructureParseRequest,
  SubmitJobAccepted,
  SubmitJobCommand,
  SupercellBuildRequest,
  SupercellBuildResponse,
  ZPEParseResponse,
  ZPEQueueTargetList,
  ZPEQueueTargetSelectResponse,
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
  exportStructureCif: (structure: Structure) => Promise<string>
  deltaTransplant: (params: DeltaTransplantRequest) => Promise<string>
  buildSupercell: (params: SupercellBuildRequest) => Promise<SupercellBuildResponse>
  parseZpeInput: (content: string, structure_id?: string | null) => Promise<ZPEParseResponse>
  fetchQueueTargets: (params?: { limit?: number; offset?: number }) => Promise<ZPEQueueTargetList>
  selectQueueTarget: (target_id: string) => Promise<ZPEQueueTargetSelectResponse>
  submitRuntimeJob: (request: SubmitJobCommand) => Promise<SubmitJobAccepted>
  postRuntimeEvent: (job_id: string, event: ExecutionEvent) => Promise<RuntimeEventAck>
  fetchRuntimeStatus: (job_id: string) => Promise<RuntimeJobStatus>
  fetchRuntimeDetail: (job_id: string) => Promise<RuntimeJobDetail>
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

    async exportStructureCif(structure) {
      const body: StructureExportRequest = { structure }
      return requestText({
        path: '/structures/cif',
        method: 'POST',
        body,
      })
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
        path: '/runtime/parse',
        method: 'POST',
        body: payload,
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
        path: `/runtime/targets${query ? `?${query}` : ''}`,
        token: withToken(),
      })
    },

    async selectQueueTarget(target_id) {
      const safeId = encodeURIComponent(target_id)
      return requestJson<ZPEQueueTargetSelectResponse>({
        path: `/runtime/targets/${safeId}/active`,
        method: 'PUT',
        token: withToken(),
      })
    },

    async submitRuntimeJob(requestPayload) {
      return requestJson<SubmitJobAccepted>({
        path: '/runtime/jobs:submit',
        method: 'POST',
        token: withToken(),
        body: requestPayload,
      })
    },

    async postRuntimeEvent(job_id, event) {
      const safeId = encodeURIComponent(job_id)
      return requestJson<RuntimeEventAck>({
        path: `/runtime/jobs/${safeId}/events`,
        method: 'POST',
        token: withToken(),
        body: event,
      })
    },

    async fetchRuntimeStatus(job_id) {
      const safeId = encodeURIComponent(job_id)
      return requestJson<RuntimeJobStatus>({
        path: `/runtime/jobs/${safeId}`,
        token: withToken(),
      })
    },

    async fetchRuntimeDetail(job_id) {
      const safeId = encodeURIComponent(job_id)
      return requestJson<RuntimeJobDetail>({
        path: `/runtime/jobs/${safeId}/detail`,
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
