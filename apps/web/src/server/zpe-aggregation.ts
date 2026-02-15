import { createServerFn } from '@tanstack/react-start'
import { JOB_STATES } from '@chem-model/shared/job-state'
import type { JobState } from '@chem-model/shared/job-state'
import type { ZPEJobStatus } from '@chem-model/api-client'

import { requestApiInternal } from '@/server/api'

export type AggregatedZpeJobStatusRequest = {
  jobId: string
  token?: string | null
}

type LegacyZpeJobStatusPayload = {
  status: JobState
  detail?: string | null
  updated_at?: string | null
}

type ConvexProjectionJobStatusPayload = {
  jobId: string
  status: JobState
  detail?: string | null
  updatedAt?: string | null
  eventTime?: string | null
}

type AdapterDetailJobStatusPayload = {
  jobId?: string
  status?: JobState
  detail?: string | null
  updated_at?: string | null
  updatedAt?: string | null
}

export type AggregatedZpeErrorCode =
  | 'AUTH_REQUIRED'
  | 'UPSTREAM_UNAVAILABLE'
  | 'UPSTREAM_PAYLOAD_UNSUPPORTED'
  | 'UPSTREAM_JOB_MISMATCH'
  | 'UPSTREAM_STATUS_MISMATCH'
  | 'UNKNOWN'

export type AggregatedZpeErrorEnvelope = {
  error: {
    code: AggregatedZpeErrorCode
    message: string
    details?: Record<string, unknown>
  }
}

export class AggregatedZpeError extends Error {
  readonly envelope: AggregatedZpeErrorEnvelope

  constructor(
    code: AggregatedZpeErrorCode,
    message: string,
    details?: Record<string, unknown>,
  ) {
    super(message)
    this.name = 'AggregatedZpeError'
    this.envelope = {
      error: {
        code,
        message,
        details,
      },
    }
  }
}

const createAggregationError = (
  code: AggregatedZpeErrorCode,
  message: string,
  details?: Record<string, unknown>,
): AggregatedZpeError => {
  return new AggregatedZpeError(code, message, details)
}

export const toAggregatedZpeErrorEnvelope = (
  error: unknown,
): AggregatedZpeErrorEnvelope => {
  if (error instanceof AggregatedZpeError) {
    return error.envelope
  }
  if (error instanceof Error) {
    return {
      error: {
        code: 'UNKNOWN',
        message: error.message,
      },
    }
  }
  return {
    error: {
      code: 'UNKNOWN',
      message: 'Unknown aggregation error',
    },
  }
}

const JOB_STATE_SET = new Set<JobState>(JOB_STATES)

const isRecord = (value: unknown): value is Record<string, unknown> => {
  return typeof value === 'object' && value !== null
}

const isStringOrNullish = (
  value: unknown,
): value is string | null | undefined => {
  return value === null || value === undefined || typeof value === 'string'
}

const isJobState = (value: unknown): value is JobState => {
  return typeof value === 'string' && JOB_STATE_SET.has(value as JobState)
}

const isLegacyStatusPayload = (
  value: unknown,
): value is LegacyZpeJobStatusPayload => {
  if (!isRecord(value)) {
    return false
  }
  return (
    isJobState(value.status) &&
    isStringOrNullish(value.detail) &&
    isStringOrNullish(value.updated_at)
  )
}

const isConvexProjectionStatusPayload = (
  value: unknown,
): value is ConvexProjectionJobStatusPayload => {
  if (!isRecord(value)) {
    return false
  }
  return (
    typeof value.jobId === 'string' &&
    isJobState(value.status) &&
    isStringOrNullish(value.detail) &&
    isStringOrNullish(value.updatedAt) &&
    isStringOrNullish(value.eventTime)
  )
}

const isAdapterDetailStatusPayload = (
  value: unknown,
): value is AdapterDetailJobStatusPayload => {
  if (!isRecord(value)) {
    return false
  }
  return (
    (value.jobId === undefined || typeof value.jobId === 'string') &&
    (value.status === undefined || isJobState(value.status)) &&
    isStringOrNullish(value.detail) &&
    isStringOrNullish(value.updated_at) &&
    isStringOrNullish(value.updatedAt)
  )
}

export const requireAuthToken = (token: string | null | undefined): string => {
  if (!token) {
    throw createAggregationError('AUTH_REQUIRED', 'Authentication required')
  }
  return token
}

export const normalizeAggregatedZpeJobStatus = (
  payload: unknown,
  expectedJobId: string,
): ZPEJobStatus => {
  const hasConvexMarkers =
    isRecord(payload) &&
    ('jobId' in payload || 'updatedAt' in payload || 'eventTime' in payload)

  if (hasConvexMarkers) {
    if (!isConvexProjectionStatusPayload(payload)) {
      throw createAggregationError(
        'UPSTREAM_PAYLOAD_UNSUPPORTED',
        'Unsupported status payload from upstream',
      )
    }
    if (payload.jobId !== expectedJobId) {
      throw createAggregationError(
        'UPSTREAM_JOB_MISMATCH',
        'Upstream job status payload mismatch',
      )
    }
    return {
      status: payload.status,
      detail: payload.detail ?? null,
      updated_at: payload.updatedAt ?? payload.eventTime ?? null,
    }
  }

  if (isLegacyStatusPayload(payload)) {
    return {
      status: payload.status,
      detail: payload.detail ?? null,
      updated_at: payload.updated_at ?? null,
    }
  }

  throw createAggregationError(
    'UPSTREAM_PAYLOAD_UNSUPPORTED',
    'Unsupported status payload from upstream',
  )
}

export const normalizeAggregatedZpeJobStatusFromSources = (
  projectionPayload: unknown,
  adapterPayload: unknown,
  expectedJobId: string,
): ZPEJobStatus => {
  const projection = normalizeAggregatedZpeJobStatus(projectionPayload, expectedJobId)

  if (!isAdapterDetailStatusPayload(adapterPayload)) {
    throw createAggregationError(
      'UPSTREAM_PAYLOAD_UNSUPPORTED',
      'Unsupported status payload from upstream',
    )
  }

  if (adapterPayload.jobId && adapterPayload.jobId !== expectedJobId) {
    throw createAggregationError(
      'UPSTREAM_JOB_MISMATCH',
      'Upstream job status payload mismatch',
    )
  }

  if (
    adapterPayload.status !== undefined &&
    adapterPayload.status !== projection.status
  ) {
    throw createAggregationError(
      'UPSTREAM_STATUS_MISMATCH',
      'Upstream status payload mismatch',
      {
        projectionStatus: projection.status,
        adapterStatus: adapterPayload.status,
      },
    )
  }

  return {
    status: projection.status,
    detail: adapterPayload.detail ?? projection.detail ?? null,
    updated_at:
      projection.updated_at ??
      adapterPayload.updated_at ??
      adapterPayload.updatedAt ??
      null,
  }
}

type UpstreamStatusRequester = (request: {
  path: string
  token: string
}) => Promise<unknown>

const DEFAULT_UPSTREAM_TIMEOUT_MS = 3000
const SECONDARY_SOURCE_GRACE_MS = 25

const withTimeout = <T>(
  promise: Promise<T>,
  timeoutMs: number,
  source: 'projection' | 'adapter',
): Promise<T> => {
  return new Promise<T>((resolve, reject) => {
    const timer = setTimeout(() => {
      reject(
        createAggregationError(
          'UPSTREAM_UNAVAILABLE',
          `Timed out fetching ${source} job status`,
        ),
      )
    }, timeoutMs)

    promise.then(
      (value) => {
        clearTimeout(timer)
        resolve(value)
      },
      (error: unknown) => {
        clearTimeout(timer)
        reject(error)
      },
    )
  })
}

export const fetchAggregatedZpeJobStatusFromUpstreams = async ({
  jobId,
  token,
  requester = requestApiInternal,
  timeoutMs = DEFAULT_UPSTREAM_TIMEOUT_MS,
}: {
  jobId: string
  token: string
  requester?: UpstreamStatusRequester
  timeoutMs?: number
}): Promise<ZPEJobStatus> => {
  const safeJobId = encodeURIComponent(jobId)
  const projectionPath = `/zpe/jobs/${safeJobId}/projection`
  const adapterStatusPath = `/zpe/jobs/${safeJobId}`

  type Source = 'projection' | 'adapter'
  type SourceResult =
    | {
        source: Source
        status: 'fulfilled'
        value: unknown
      }
    | {
        source: Source
        status: 'rejected'
        reason: unknown
      }
  type PendingResult = {
    status: 'pending'
  }

  const projectionPromise: Promise<SourceResult> = withTimeout(
    requester({
      path: projectionPath,
      token,
    }),
    timeoutMs,
    'projection',
  ).then(
    (value) => ({
      source: 'projection',
      status: 'fulfilled',
      value,
    }),
    (reason: unknown) => ({
      source: 'projection',
      status: 'rejected',
      reason,
    }),
  )

  const adapterPromise: Promise<SourceResult> = withTimeout(
    requester({
      path: adapterStatusPath,
      token,
    }),
    timeoutMs,
    'adapter',
  ).then(
    (value) => ({
      source: 'adapter',
      status: 'fulfilled',
      value,
    }),
    (reason: unknown) => ({
      source: 'adapter',
      status: 'rejected',
      reason,
    }),
  )

  const firstSettled = await Promise.race([projectionPromise, adapterPromise])
  const secondaryPromise =
    firstSettled.source === 'projection' ? adapterPromise : projectionPromise

  const secondarySettled = await Promise.race<SourceResult | PendingResult>([
    secondaryPromise,
    new Promise<PendingResult>((resolve) => {
      setTimeout(() => resolve({ status: 'pending' }), SECONDARY_SOURCE_GRACE_MS)
    }),
  ])

  const settled: Array<SourceResult> = [
    firstSettled,
    ...(secondarySettled.status === 'pending' ? [] : [secondarySettled]),
  ]

  const projectionResult = settled.find(
    (result) => result.source === 'projection',
  )
  const adapterResult = settled.find((result) => result.source === 'adapter')

  if (
    projectionResult?.status === 'fulfilled' &&
    adapterResult?.status === 'fulfilled'
  ) {
    return normalizeAggregatedZpeJobStatusFromSources(
      projectionResult.value,
      adapterResult.value,
      jobId,
    )
  }

  if (projectionResult?.status === 'fulfilled') {
    return normalizeAggregatedZpeJobStatus(projectionResult.value, jobId)
  }

  if (adapterResult?.status === 'fulfilled') {
    return normalizeAggregatedZpeJobStatus(adapterResult.value, jobId)
  }

  if (secondarySettled.status === 'pending') {
    const finalResult = await secondaryPromise
    if (finalResult.status === 'fulfilled') {
      return normalizeAggregatedZpeJobStatus(finalResult.value, jobId)
    }
  }

  throw createAggregationError(
    'UPSTREAM_UNAVAILABLE',
    'Failed to fetch upstream job status',
  )
}

type FetchAggregatedZpeJobStatusFn = ((options: {
  data: AggregatedZpeJobStatusRequest
}) => Promise<ZPEJobStatus>) & {
  url: string
}

export const fetchAggregatedZpeJobStatus: FetchAggregatedZpeJobStatusFn =
  createServerFn({
    method: 'POST',
  })
    .inputValidator((data: AggregatedZpeJobStatusRequest) => data)
    .handler(async ({ data }) => {
      const token = requireAuthToken(data.token)
      return fetchAggregatedZpeJobStatusFromUpstreams({
        jobId: data.jobId,
        token,
      })
    })
