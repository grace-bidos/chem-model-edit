import { createServerFn } from '@tanstack/react-start'
import { JOB_STATES } from '@chem-model/shared/job-state'
import type { JobState } from '@chem-model/shared/job-state'
import type { ZPEJobStatus } from '@chem-model/api-client'

import { requestApiInternal } from '@/server/api'

export type AggregatedZpeJobStatusRequest = {
  jobId: string
  token?: string | null
}

type ProjectionJobStatusPayload = {
  status: JobState
  detail?: string | null
  updated_at?: string | null
}

export type AggregatedZpeErrorCode =
  | 'AUTH_REQUIRED'
  | 'UPSTREAM_UNAVAILABLE'
  | 'UPSTREAM_JOB_MISMATCH'
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

const isJobState = (value: unknown): value is JobState => {
  return typeof value === 'string' && JOB_STATE_SET.has(value as JobState)
}

const isProjectionStatusPayload = (
  value: unknown,
): value is ProjectionJobStatusPayload => {
  if (!isRecord(value)) {
    return false
  }
  return (
    isJobState(value.status) &&
    (value.detail === undefined || value.detail === null || typeof value.detail === 'string') &&
    (value.updated_at === undefined ||
      value.updated_at === null ||
      typeof value.updated_at === 'string')
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
): ZPEJobStatus => {
  if (!isProjectionStatusPayload(payload)) {
    throw createAggregationError(
      'UPSTREAM_UNAVAILABLE',
      'Unsupported projection payload from upstream',
    )
  }

  return {
    status: payload.status,
    detail: payload.detail ?? null,
    updated_at: payload.updated_at ?? null,
  }
}

type UpstreamStatusRequester = (request: {
  path: string
  token: string
}) => Promise<unknown>

const DEFAULT_UPSTREAM_TIMEOUT_MS = 3000

const withTimeout = <T>(
  promise: Promise<T>,
  timeoutMs: number,
): Promise<T> => {
  return new Promise<T>((resolve, reject) => {
    const timer = setTimeout(() => {
      reject(
        createAggregationError(
          'UPSTREAM_UNAVAILABLE',
          'Timed out fetching projection job status',
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
  const projectionPath = `/runtime/jobs/${safeJobId}/projection`

  try {
    const payload = await withTimeout(
      requester({
        path: projectionPath,
        token,
      }),
      timeoutMs,
    )
    return normalizeAggregatedZpeJobStatus(payload)
  } catch (error: unknown) {
    if (error instanceof AggregatedZpeError) {
      throw error
    }
    throw createAggregationError(
      'UPSTREAM_UNAVAILABLE',
      'Failed to fetch projection job status',
      {
        jobId,
      },
    )
  }
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
