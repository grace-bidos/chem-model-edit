import { describe, expect, it, vi } from 'vitest'

import {
  AggregatedZpeError,
  fetchAggregatedZpeJobStatusFromUpstreams,
  normalizeAggregatedZpeJobStatus,
  requireAuthToken,
  toAggregatedZpeErrorEnvelope,
} from './zpe-aggregation'

describe('requireAuthToken', () => {
  it('throws when token is missing', () => {
    expect(() => requireAuthToken(null)).toThrow('Authentication required')
  })

  it('returns token when present', () => {
    expect(requireAuthToken('token-123')).toBe('token-123')
  })
})

describe('normalizeAggregatedZpeJobStatus', () => {
  it('accepts runtime projection payloads', () => {
    const normalized = normalizeAggregatedZpeJobStatus({
      status: 'started',
      detail: 'running',
      updated_at: '2026-02-15T08:00:00Z',
    })

    expect(normalized).toEqual({
      status: 'started',
      detail: 'running',
      updated_at: '2026-02-15T08:00:00Z',
    })
  })

  it('rejects unsupported payloads', () => {
    expect(() =>
      normalizeAggregatedZpeJobStatus({
        detail: 'missing-status',
      }),
    ).toThrow('Unsupported projection payload from upstream')
  })
})

describe('fetchAggregatedZpeJobStatusFromUpstreams', () => {
  it('uses projection endpoint only', async () => {
    const calls: Array<string> = []
    const requester = (request: { path: string; token: string }) => {
      calls.push(`${request.path}|${request.token}`)
      return Promise.resolve({
        status: 'finished',
        detail: 'done',
        updated_at: '2026-02-15T08:20:00Z',
      })
    }

    const normalized = await fetchAggregatedZpeJobStatusFromUpstreams({
      jobId: 'job-9',
      token: 'token-9',
      requester,
    })

    expect(calls).toEqual(['/runtime/jobs/job-9/projection|token-9'])
    expect(normalized).toEqual({
      status: 'finished',
      detail: 'done',
      updated_at: '2026-02-15T08:20:00Z',
    })
  })

  it('times out when projection source hangs', async () => {
    vi.useFakeTimers()
    try {
      const request = fetchAggregatedZpeJobStatusFromUpstreams({
        jobId: 'job-13',
        token: 'token-13',
        timeoutMs: 25,
        requester: () => new Promise<unknown>(() => {}),
      })
      const assertion = expect(request).rejects.toThrow(
        'Timed out fetching projection job status',
      )

      await vi.advanceTimersByTimeAsync(25)
      await assertion
    } finally {
      vi.useRealTimers()
    }
  })

  it('throws unavailable when requester fails', async () => {
    await expect(
      fetchAggregatedZpeJobStatusFromUpstreams({
        jobId: 'job-12',
        token: 'token-12',
        requester: () => Promise.reject(new Error('down')),
      }),
    ).rejects.toThrow('Failed to fetch projection job status')
  })
})

describe('toAggregatedZpeErrorEnvelope', () => {
  it('serializes AggregatedZpeError', () => {
    const error = new AggregatedZpeError(
      'UPSTREAM_UNAVAILABLE',
      'Projection unavailable',
      { jobId: 'job-1' },
    )
    expect(toAggregatedZpeErrorEnvelope(error)).toEqual({
      error: {
        code: 'UPSTREAM_UNAVAILABLE',
        message: 'Projection unavailable',
        details: { jobId: 'job-1' },
      },
    })
  })

  it('serializes native Error as UNKNOWN', () => {
    const envelope = toAggregatedZpeErrorEnvelope(new Error('boom'))
    expect(envelope.error.code).toBe('UNKNOWN')
    expect(envelope.error.message).toBe('boom')
  })
})
