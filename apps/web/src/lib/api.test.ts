/* @vitest-environment jsdom */
import { beforeEach, describe, expect, it, vi } from 'vitest'

import {
  fetchZpeStatus,
  parseQeInput,
  setApiTokenProvider,
  structureViewUrl,
} from './api'
import { fetchAggregatedZpeJobStatus } from '@/server/zpe-aggregation'
import { requestApi } from '@/server/api'

vi.mock('@/server/api', () => ({
  requestApi: vi.fn(),
}))

vi.mock('@/server/zpe-aggregation', () => ({
  fetchAggregatedZpeJobStatus: vi.fn(),
}))

vi.mock('@chem-model/api-client', () => ({
  createApiClient: ({ request }: { request: (params: unknown) => Promise<unknown> }) => ({
    parseQeInput: (body: unknown) =>
      request({
        path: '/qe/parse',
        method: 'POST',
        body,
      }),
    createStructureFromQe: vi.fn(),
    getStructure: vi.fn(),
    exportQeInput: vi.fn(),
    exportStructureCif: vi.fn(),
    deltaTransplant: vi.fn(),
    buildSupercell: vi.fn(),
    parseZpeInput: vi.fn(),
    createEnrollToken: vi.fn(),
    fetchQueueTargets: vi.fn(),
    selectQueueTarget: vi.fn(),
    createZpeJob: vi.fn(),
    fetchZpeResult: vi.fn(),
    downloadZpeFile: vi.fn(),
    structureViewPath: (structureId: string, params?: { format?: 'cif' }) => {
      const query = params?.format ? `?format=${params.format}` : ''
      return `/structures/${structureId}${query}`
    },
  }),
}))

const mockedRequestApi = vi.mocked(requestApi)
const mockedFetchAggregated = vi.mocked(fetchAggregatedZpeJobStatus)

describe('api client integration helpers', () => {
  beforeEach(() => {
    vi.clearAllMocks()
    setApiTokenProvider(null)
    vi.unstubAllEnvs()
    delete window.__API_BASE__
  })

  it('uses token provider when request token is not specified', async () => {
    mockedRequestApi.mockResolvedValueOnce({ ok: true })
    setApiTokenProvider(() => Promise.resolve('from-provider'))

    await parseQeInput('&CONTROL\n/')

    expect(mockedRequestApi).toHaveBeenCalledWith({
      data: {
        path: '/qe/parse',
        method: 'POST',
        body: '&CONTROL\n/',
        token: 'from-provider',
      },
    })
  })

  it('passes null token when fetching zpe status without provider token', async () => {
    mockedFetchAggregated.mockResolvedValueOnce({
      status: 'queued',
      detail: null,
      updated_at: null,
    })

    await fetchZpeStatus('job-42')

    expect(mockedFetchAggregated).toHaveBeenCalledWith({
      data: {
        jobId: 'job-42',
        token: null,
      },
    })
  })

  it('uses provider token when fetching zpe status', async () => {
    mockedFetchAggregated.mockResolvedValueOnce({
      status: 'started',
      detail: 'running',
      updated_at: null,
    })
    setApiTokenProvider(() => Promise.resolve('status-token'))

    await fetchZpeStatus('job-99')

    expect(mockedFetchAggregated).toHaveBeenCalledWith({
      data: {
        jobId: 'job-99',
        token: 'status-token',
      },
    })
  })

  it('normalizes api base from window.__API_BASE__ and appends structure path', () => {
    window.__API_BASE__ = 'https://example.test/root/'

    expect(structureViewUrl('abc')).toBe('https://example.test/root/api/structures/abc')
  })

  it('keeps /api suffix and includes query params for structure view url', () => {
    window.__API_BASE__ = 'https://example.test/custom/api'

    expect(structureViewUrl('def', { format: 'cif' })).toBe(
      'https://example.test/custom/api/structures/def?format=cif',
    )
  })

  it('falls back to VITE_API_BASE when window override is missing', () => {
    vi.stubEnv('VITE_API_BASE', 'https://api.example.test/base/')

    expect(structureViewUrl('ghi')).toBe('https://api.example.test/base/api/structures/ghi')
  })
})
