import { describe, expect, it } from 'vitest'

import {
  normalizeAggregatedZpeJobStatus,
  normalizeAggregatedZpeJobStatusFromSources,
  toAggregatedZpeErrorEnvelope,
} from './zpe-aggregation'

type ContractFixture = {
  name: string
  expectedJobId: string
  projection: unknown
  adapter?: unknown
  expected: {
    status: 'queued' | 'started' | 'finished' | 'failed'
    detail: string | null
    updated_at: string | null
  }
}

const CONTRACT_FIXTURES: ContractFixture[] = [
  {
    name: 'legacy FastAPI payload remains compatible',
    expectedJobId: 'job-contract-1',
    projection: {
      status: 'started',
      detail: 'running',
      updated_at: '2026-02-15T09:00:00Z',
    },
    expected: {
      status: 'started',
      detail: 'running',
      updated_at: '2026-02-15T09:00:00Z',
    },
  },
  {
    name: 'convex projection payload maps camelCase timestamp',
    expectedJobId: 'job-contract-2',
    projection: {
      jobId: 'job-contract-2',
      status: 'finished',
      detail: null,
      updatedAt: '2026-02-15T09:01:00Z',
    },
    expected: {
      status: 'finished',
      detail: null,
      updated_at: '2026-02-15T09:01:00Z',
    },
  },
  {
    name: 'projection plus adapter detail prefers adapter detail',
    expectedJobId: 'job-contract-3',
    projection: {
      jobId: 'job-contract-3',
      status: 'queued',
      eventTime: '2026-02-15T09:02:00Z',
    },
    adapter: {
      jobId: 'job-contract-3',
      status: 'queued',
      detail: 'queued on adapter',
      updated_at: '2026-02-15T09:01:59Z',
    },
    expected: {
      status: 'queued',
      detail: 'queued on adapter',
      updated_at: '2026-02-15T09:02:00Z',
    },
  },
]

describe('Convex/BFF contract fixtures', () => {
  for (const fixture of CONTRACT_FIXTURES) {
    it(fixture.name, () => {
      const normalized = fixture.adapter
        ? normalizeAggregatedZpeJobStatusFromSources(
            fixture.projection,
            fixture.adapter,
            fixture.expectedJobId,
          )
        : normalizeAggregatedZpeJobStatus(
            fixture.projection,
            fixture.expectedJobId,
          )

      expect(normalized).toEqual(fixture.expected)
    })
  }

  it('returns compatible typed envelope for projection/adapter mismatch errors', () => {
    const envelope = toAggregatedZpeErrorEnvelope(
      (() => {
        try {
          normalizeAggregatedZpeJobStatusFromSources(
            {
              jobId: 'job-contract-4',
              status: 'started',
              updatedAt: '2026-02-15T09:03:00Z',
            },
            {
              jobId: 'job-contract-4',
              status: 'failed',
              detail: 'worker failure',
            },
            'job-contract-4',
          )
          throw new Error('expected mismatch')
        } catch (error) {
          return error
        }
      })(),
    )

    expect(envelope.error.code).toBe('UPSTREAM_STATUS_MISMATCH')
    expect(envelope.error.message).toBe('Upstream status payload mismatch')
  })
})
