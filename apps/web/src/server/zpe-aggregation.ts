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

const JOB_STATE_SET = new Set<JobState>(JOB_STATES)

const isRecord = (value: unknown): value is Record<string, unknown> => {
  return typeof value === 'object' && value !== null
}

const isStringOrNullish = (value: unknown): value is string | null | undefined => {
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

export const requireAuthToken = (token: string | null | undefined): string => {
  if (!token) {
    throw new Error('Authentication required')
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
      throw new Error('Unsupported status payload from upstream')
    }
    if (payload.jobId !== expectedJobId) {
      throw new Error('Upstream job status payload mismatch')
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

  throw new Error('Unsupported status payload from upstream')
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
      const safeJobId = encodeURIComponent(data.jobId)

      // TODO(GRA-18): fan-out to Convex projection and AiiDA detail sources.
      const upstreamPayload = await requestApiInternal<unknown>({
        path: `/zpe/jobs/${safeJobId}`,
        token,
      })

      return normalizeAggregatedZpeJobStatus(upstreamPayload, data.jobId)
    })
