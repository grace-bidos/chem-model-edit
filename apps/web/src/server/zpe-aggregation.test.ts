import { describe, expect, it } from 'vitest'

import {
  normalizeAggregatedZpeJobStatus,
  requireAuthToken,
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
