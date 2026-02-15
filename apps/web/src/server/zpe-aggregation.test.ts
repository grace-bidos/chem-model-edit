import { describe, expect, it } from 'vitest'

import {
  AggregatedZpeError,
  fetchAggregatedZpeJobStatusFromUpstreams,
  normalizeAggregatedZpeJobStatus,
  normalizeAggregatedZpeJobStatusFromSources,
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
  it('normalizes legacy FastAPI status payloads', () => {
    const normalized = normalizeAggregatedZpeJobStatus(
      {
        status: 'started',
        detail: 'running',
        updated_at: '2026-02-15T08:00:00Z',
      },
      'job-1',
    )

    expect(normalized).toEqual({
      status: 'started',
      detail: 'running',
      updated_at: '2026-02-15T08:00:00Z',
    })
  })

  it('normalizes Convex-style status payloads', () => {
    const normalized = normalizeAggregatedZpeJobStatus(
      {
        jobId: 'job-2',
        status: 'finished',
        detail: null,
        updatedAt: '2026-02-15T08:01:00Z',
      },
      'job-2',
    )

    expect(normalized).toEqual({
      status: 'finished',
      detail: null,
      updated_at: '2026-02-15T08:01:00Z',
    })
  })

  it('falls back to eventTime when updatedAt is missing', () => {
    const normalized = normalizeAggregatedZpeJobStatus(
      {
        jobId: 'job-3',
        status: 'queued',
        eventTime: '2026-02-15T08:02:00Z',
      },
      'job-3',
    )

    expect(normalized).toEqual({
      status: 'queued',
      detail: null,
      updated_at: '2026-02-15T08:02:00Z',
    })
  })

  it('rejects job id mismatch for Convex-style payloads', () => {
    expect(() =>
      normalizeAggregatedZpeJobStatus(
        {
          jobId: 'other-job',
          status: 'queued',
        },
        'job-4',
      ),
    ).toThrow('Upstream job status payload mismatch')
  })

  it('rejects unsupported payloads', () => {
    expect(() =>
      normalizeAggregatedZpeJobStatus(
        {
          status: 'queued',
          updatedAt: '2026-02-15T08:03:00Z',
        },
        'job-5',
      ),
    ).toThrow('Unsupported status payload from upstream')
  })
})

describe('normalizeAggregatedZpeJobStatusFromSources', () => {
  it('merges projection status with adapter detail', () => {
    const normalized = normalizeAggregatedZpeJobStatusFromSources(
      {
        jobId: 'job-6',
        status: 'started',
        updatedAt: '2026-02-15T08:10:00Z',
      },
      {
        jobId: 'job-6',
        status: 'started',
        detail: 'adapter-running',
        updated_at: '2026-02-15T08:09:59Z',
      },
      'job-6',
    )

    expect(normalized).toEqual({
      status: 'started',
      detail: 'adapter-running',
      updated_at: '2026-02-15T08:10:00Z',
    })
  })

  it('keeps compatibility by falling back to projection detail when adapter detail missing', () => {
    const normalized = normalizeAggregatedZpeJobStatusFromSources(
      {
        jobId: 'job-7',
        status: 'queued',
        detail: 'projection-queued',
        eventTime: '2026-02-15T08:11:00Z',
      },
      {
        jobId: 'job-7',
        status: 'queued',
      },
      'job-7',
    )

    expect(normalized).toEqual({
      status: 'queued',
      detail: 'projection-queued',
      updated_at: '2026-02-15T08:11:00Z',
    })
  })

  it('rejects status mismatch across projection and adapter', () => {
    expect(() =>
      normalizeAggregatedZpeJobStatusFromSources(
        {
          jobId: 'job-8',
          status: 'started',
          updatedAt: '2026-02-15T08:12:00Z',
        },
        {
          jobId: 'job-8',
          status: 'failed',
          detail: 'boom',
        },
        'job-8',
      ),
    ).toThrow('Upstream status payload mismatch')
  })
})

describe('fetchAggregatedZpeJobStatusFromUpstreams', () => {
  it('fans out to projection and adapter endpoints, then merges', async () => {
    const calls: Array<string> = []
    const requester = (request: { path: string; token: string }) => {
      calls.push(`${request.path}|${request.token}`)
      if (request.path.endsWith('/projection')) {
        return Promise.resolve({
          jobId: 'job-9',
          status: 'finished',
          updatedAt: '2026-02-15T08:20:00Z',
        })
      }
      return Promise.resolve({
        jobId: 'job-9',
        status: 'finished',
        detail: 'done',
        updated_at: '2026-02-15T08:19:59Z',
      })
    }

    const normalized = await fetchAggregatedZpeJobStatusFromUpstreams({
      jobId: 'job-9',
      token: 'token-9',
      requester,
    })

    expect(calls).toEqual([
      '/zpe/jobs/job-9/projection|token-9',
      '/zpe/jobs/job-9|token-9',
    ])
    expect(normalized).toEqual({
      status: 'finished',
      detail: 'done',
      updated_at: '2026-02-15T08:20:00Z',
    })
  })

  it('falls back to adapter payload when projection source fails', async () => {
    const requester = (request: { path: string }) => {
      if (request.path.endsWith('/projection')) {
        return Promise.reject(new Error('projection unavailable'))
      }
      return Promise.resolve({
        status: 'started',
        detail: 'legacy-running',
        updated_at: '2026-02-15T08:21:00Z',
      })
    }

    const normalized = await fetchAggregatedZpeJobStatusFromUpstreams({
      jobId: 'job-10',
      token: 'token-10',
      requester,
    })

    expect(normalized).toEqual({
      status: 'started',
      detail: 'legacy-running',
      updated_at: '2026-02-15T08:21:00Z',
    })
  })

  it('falls back to projection payload when adapter source fails', async () => {
    const requester = (request: { path: string }) => {
      if (!request.path.endsWith('/projection')) {
        return Promise.reject(new Error('adapter unavailable'))
      }
      return Promise.resolve({
        jobId: 'job-11',
        status: 'queued',
        eventTime: '2026-02-15T08:22:00Z',
      })
    }

    const normalized = await fetchAggregatedZpeJobStatusFromUpstreams({
      jobId: 'job-11',
      token: 'token-11',
      requester,
    })

    expect(normalized).toEqual({
      status: 'queued',
      detail: null,
      updated_at: '2026-02-15T08:22:00Z',
    })
  })

  it('throws typed unavailable error when both sources fail', async () => {
    await expect(
      fetchAggregatedZpeJobStatusFromUpstreams({
        jobId: 'job-12',
        token: 'token-12',
        requester: () => Promise.reject(new Error('down')),
      }),
    ).rejects.toThrow('Failed to fetch upstream job status')
  })
})

describe('toAggregatedZpeErrorEnvelope', () => {
  it('returns envelope from AggregatedZpeError', () => {
    const envelope = toAggregatedZpeErrorEnvelope(
      new AggregatedZpeError('UPSTREAM_JOB_MISMATCH', 'payload mismatch', {
        expectedJobId: 'job-1',
      }),
    )

    expect(envelope).toEqual({
      error: {
        code: 'UPSTREAM_JOB_MISMATCH',
        message: 'payload mismatch',
        details: {
          expectedJobId: 'job-1',
        },
      },
    })
  })

  it('maps unknown errors to UNKNOWN envelope shape', () => {
    const envelope = toAggregatedZpeErrorEnvelope(new Error('boom'))
    expect(envelope).toEqual({
      error: {
        code: 'UNKNOWN',
        message: 'boom',
      },
    })
  })
})
