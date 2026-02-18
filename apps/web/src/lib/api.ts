import { createApiClient } from '@chem-model/api-client'

import type {
  ApiRequest,
  SubmitJobCommand,
  ZPEJobStatus,
} from '@chem-model/api-client'
import { requestApi } from '@/server/api'
import { fetchAggregatedZpeJobStatus } from '@/server/zpe-aggregation'

type TokenProvider = () => Promise<string | null>

let tokenProvider: TokenProvider | null = null

export function setApiTokenProvider(provider: TokenProvider | null) {
  tokenProvider = provider
}

const request = async <T>(params: ApiRequest): Promise<T> => {
  const token =
    params.token === undefined
      ? ((await tokenProvider?.()) ?? undefined)
      : params.token
  return (await requestApi({ data: { ...params, token } })) as T
}

const api = createApiClient({
  request,
})

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

const decodeJwtPayload = (token: string): Record<string, unknown> | null => {
  const parts = token.split('.')
  if (parts.length < 2) {
    return null
  }
  try {
    const payloadPart = parts[1] ?? ''
    const normalized = payloadPart.replace(/-/g, '+').replace(/_/g, '/')
    const padded = normalized.padEnd(
      normalized.length + ((4 - (normalized.length % 4)) % 4),
      '=',
    )
    const payload = JSON.parse(atob(padded)) as unknown
    return payload && typeof payload === 'object'
      ? (payload as Record<string, unknown>)
      : null
  } catch (_error) {
    return null
  }
}

const resolveClaimString = (
  claims: Record<string, unknown> | null,
  keys: Array<string>,
): string | null => {
  if (!claims) {
    return null
  }
  for (const key of keys) {
    const value = claims[key]
    if (typeof value === 'string' && value.trim().length > 0) {
      return value.trim()
    }
  }
  return null
}

export type RuntimeSubmitInput = {
  queueName: string
  managementNodeId: string
}

export const parseQeInput = api.parseQeInput
export const createStructureFromQe = api.createStructureFromQe
export const getStructure = api.getStructure
export const exportQeInput = api.exportQeInput
export const exportStructureCif = api.exportStructureCif
export const deltaTransplant = api.deltaTransplant
export const buildSupercell = api.buildSupercell
export const parseZpeInput = api.parseZpeInput
export const fetchQueueTargets = api.fetchQueueTargets
export const selectQueueTarget = api.selectQueueTarget

export const createZpeJob = async (
  input: RuntimeSubmitInput,
): Promise<{ id: string }> => {
  const token = (await tokenProvider?.()) ?? null
  if (!token) {
    throw new Error('Authentication required')
  }
  const claims = decodeJwtPayload(token)
  const tenantId =
    resolveClaimString(claims, ['tenant_id', 'tenantId', 'org_id', 'orgId']) ??
    null
  const userId = resolveClaimString(claims, ['sub', 'user_id', 'userId']) ?? null
  if (!tenantId || !userId) {
    throw new Error('Missing required tenant/user claims for runtime submit')
  }

  const jobId = `job-${crypto.randomUUID()}`
  const command: SubmitJobCommand = {
    tenant_id: tenantId,
    workspace_id: tenantId,
    job_id: jobId,
    idempotency_key: `submit-${jobId}`,
    management_node_id: input.managementNodeId,
    execution_profile: {
      queue_name: input.queueName,
      qos: null,
      account: null,
    },
    resource_shape: {
      cpu: 1,
      memory_mib: 2048,
      walltime_seconds: 3600,
    },
    payload_ref: {
      input_uri: `inline://runtime-submit/${jobId}`,
      artifact_bucket: null,
    },
    requested_by: {
      user_id: userId,
      session_id: null,
    },
  }

  const accepted = await api.submitRuntimeJob(command)
  return { id: accepted.job_id }
}

export const fetchZpeStatus = async (job_id: string): Promise<ZPEJobStatus> => {
  const token = (await tokenProvider?.()) ?? null
  return fetchAggregatedZpeJobStatus({
    data: {
      jobId: job_id,
      token,
    },
  })
}

export function structureViewUrl(
  structureId: string,
  params?: {
    format?: 'cif'
  },
): string {
  return `${resolveApiBase()}${api.structureViewPath(structureId, params)}`
}
