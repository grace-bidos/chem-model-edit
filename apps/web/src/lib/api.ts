import { createApiClient } from '@chem-model/api-client'

import { getAuthToken } from './auth'

import type { ApiRequest } from '@chem-model/api-client'
import { requestApi } from '@/server/api'

const request = async <T>(params: ApiRequest): Promise<T> => {
  return (await requestApi({ data: params })) as T
}

const api = createApiClient({
  request,
  getToken: () => getAuthToken() ?? undefined,
})

// API_BASE should be a host root or already end with "/api" (no "/api/v1" style path).
const normalizeApiBase = (base: string) => {
  const trimmed = base.endsWith('/') ? base.slice(0, -1) : base
  return trimmed.endsWith('/api') ? trimmed : `${trimmed}/api`
}

const resolveApiBase = (): string => {
  if (typeof window !== 'undefined' && window.__API_BASE__) {
    return normalizeApiBase(window.__API_BASE__)
  }
  return normalizeApiBase(
    import.meta.env.VITE_API_BASE ?? 'http://localhost:8000',
  )
}

export const parseQeInput = api.parseQeInput
export const createStructureFromQe = api.createStructureFromQe
export const getStructure = api.getStructure
export const exportQeInput = api.exportQeInput
export const exportStructureCif = api.exportStructureCif
export const deltaTransplant = api.deltaTransplant
export const buildSupercell = api.buildSupercell
export const parseZpeInput = api.parseZpeInput
export const registerAccount = api.registerAccount
export const loginAccount = api.loginAccount
export const fetchAuthMe = api.fetchAuthMe
export const logoutAccount = api.logoutAccount
export const createEnrollToken = api.createEnrollToken
export const fetchQueueTargets = api.fetchQueueTargets
export const selectQueueTarget = api.selectQueueTarget
export const createZpeJob = api.createZpeJob
export const fetchZpeStatus = api.fetchZpeStatus
export const fetchZpeResult = api.fetchZpeResult
export const downloadZpeFile = api.downloadZpeFile

export function structureViewUrl(
  structureId: string,
  params?: {
    format?: 'cif'
  },
): string {
  return `${resolveApiBase()}${api.structureViewPath(structureId, params)}`
}
