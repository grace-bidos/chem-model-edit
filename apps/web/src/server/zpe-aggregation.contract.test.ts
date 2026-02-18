import { describe, expect, it } from 'vitest'

import {
  normalizeAggregatedZpeJobStatus,
  toAggregatedZpeErrorEnvelope,
} from './zpe-aggregation'

type ContractFixture = {
  name: string
  projection: unknown
  expected: {
    status: 'queued' | 'started' | 'finished' | 'failed'
    detail: string | null
    updated_at: string | null
  }
}

const CONTRACT_FIXTURES: Array<ContractFixture> = [
  {
    name: 'runtime projection payload remains stable',
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
    name: 'runtime projection payload allows null detail',
    projection: {
      status: 'finished',
      detail: null,
      updated_at: '2026-02-15T09:01:00Z',
    },
    expected: {
      status: 'finished',
      detail: null,
      updated_at: '2026-02-15T09:01:00Z',
    },
  },
]

describe('Projection/BFF contract fixtures', () => {
  for (const fixture of CONTRACT_FIXTURES) {
    it(fixture.name, () => {
      expect(normalizeAggregatedZpeJobStatus(fixture.projection)).toEqual(
        fixture.expected,
      )
    })
  }

  it('returns typed envelope for unsupported projection payload', () => {
    const envelope = toAggregatedZpeErrorEnvelope(
      (() => {
        try {
          normalizeAggregatedZpeJobStatus({ detail: 'missing-status' })
          throw new Error('expected unsupported payload')
        } catch (error) {
          return error
        }
      })(),
    )

    expect(envelope.error.code).toBe('UPSTREAM_UNAVAILABLE')
    expect(envelope.error.message).toBe('Unsupported projection payload from upstream')
  })
})
