import { createApiClient } from '@chem-model/api-client'

import { getAuthToken } from './auth'

import type { ApiRequest } from '@chem-model/api-client'
import { requestApi } from '@/server/api'

type ApiError = {
  error?: {
    code?: string
    message?: string
    details?: unknown
  }
}

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

const requestViaFetch = async <T>(params: ApiRequest): Promise<T> => {
  const { path, method, body, token, responseType } = params
  const resolvedMethod = method ?? 'GET'
  const headers: Record<string, string> = {}
  if (body !== undefined) {
    headers['Content-Type'] = 'application/json'
  }
  if (token) {
    headers.Authorization = `Bearer ${token}`
  }
  const response = await fetch(`${resolveApiBase()}${path}`, {
    method: resolvedMethod,
    headers,
    body: body !== undefined ? JSON.stringify(body) : undefined,
  })
  if (!response.ok) {
    let message = response.statusText
    try {
      const data = (await response.clone().json()) as ApiError
      if (data.error?.message) {
        message = data.error.message
      }
    } catch (_err) {
      // ignore JSON parsing errors
    }
    throw new Error(message)
  }
  if (responseType === 'text') {
    return (await response.text()) as T
  }
  return (await response.json()) as T
}

const request = async <T>(params: ApiRequest): Promise<T> => {
  if (typeof window !== 'undefined') {
    return requestViaFetch<T>(params)
  }
  return (await requestApi({ data: params })) as T
}

const api = createApiClient({
  request,
  getToken: () => getAuthToken() ?? undefined,
})

export const parseQeInput = api.parseQeInput
export const createStructureFromQe = api.createStructureFromQe
export const getStructure = api.getStructure
export const exportQeInput = api.exportQeInput
export const deltaTransplant = api.deltaTransplant
export const generateSupercell = api.generateSupercell
export const generateTiledSupercell = api.generateTiledSupercell
export const buildSupercell = api.buildSupercell
export const latticeVectorsToParams = api.latticeVectorsToParams
export const latticeParamsToVectors = api.latticeParamsToVectors
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
  const format = params?.format ?? 'cif'
  const query = new URLSearchParams({
    format,
  })
  const safeId = encodeURIComponent(structureId)
  return `${resolveApiBase()}/structures/${safeId}/view?${query.toString()}`
}
