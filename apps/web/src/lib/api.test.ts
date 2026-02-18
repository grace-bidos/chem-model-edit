/* @vitest-environment jsdom */
import { beforeEach, describe, expect, it, vi } from 'vitest'

import {
  createZpeJob,
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
    parseQeInput: (content: string) =>
      request({
        path: '/structures/parse',
        method: 'POST',
        body: { content },
      }),
    createStructureFromQe: vi.fn(),
    getStructure: vi.fn(),
    exportQeInput: vi.fn(),
    exportStructureCif: vi.fn(),
    deltaTransplant: vi.fn(),
    buildSupercell: vi.fn(),
    parseZpeInput: vi.fn(),
    fetchQueueTargets: vi.fn(),
    selectQueueTarget: vi.fn(),
    submitRuntimeJob: (body: unknown) =>
      request({
        path: '/runtime/jobs:submit',
        method: 'POST',
        body,
      }),
    postRuntimeEvent: vi.fn(),
    fetchRuntimeStatus: vi.fn(),
    fetchRuntimeDetail: vi.fn(),
    structureViewPath: (structureId: string, params?: { format?: 'cif' }) => {
      const query = params?.format ? `?format=${params.format}` : ''
      return `/structures/${structureId}${query}`
    },
  }),
}))

const mockedRequestApi = vi.mocked(requestApi)
const mockedFetchAggregated = vi.mocked(fetchAggregatedZpeJobStatus)
const makeJwt = (payload: Record<string, unknown>) => {
  const encoded = btoa(JSON.stringify(payload))
    .replace(/\+/g, '-')
    .replace(/\//g, '_')
    .replace(/=+$/g, '')
  return `header.${encoded}.sig`
}

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
        path: '/structures/parse',
        method: 'POST',
        body: { content: '&CONTROL\n/' },
        token: 'from-provider',
      },
    })
  })

  it('builds runtime submit command from JWT claims', async () => {
    mockedRequestApi.mockResolvedValueOnce({
      job_id: 'job-accepted',
      submission_id: 'sub-1',
      execution_owner: 'aiida-slurm',
      accepted_at: '2026-02-18T00:00:00Z',
      trace_id: 'trace-1',
    })
    setApiTokenProvider(() =>
      Promise.resolve(
        makeJwt({
          tenant_id: 'tenant-a',
          sub: 'user-a',
        }),
      ),
    )

    const result = await createZpeJob({
      queueName: 'default',
      managementNodeId: 'mgmt-node-1',
    })

    expect(result.id).toBe('job-accepted')
    const submitCall = mockedRequestApi.mock.calls.at(-1)?.[0] as {
      data: { path: string; method: string; body: Record<string, unknown>; token: string }
    }
    expect(submitCall.data.path).toBe('/runtime/jobs:submit')
    expect(submitCall.data.method).toBe('POST')
    expect(submitCall.data.body.tenant_id).toBe('tenant-a')
    expect(submitCall.data.body.workspace_id).toBe('tenant-a')
    expect(submitCall.data.body.management_node_id).toBe('mgmt-node-1')
    expect(submitCall.data.body.requested_by).toMatchObject({ user_id: 'user-a' })
    expect((submitCall.data.body.execution_profile as Record<string, unknown>).queue_name).toBe('default')
  })

  it('throws when required claims are missing for runtime submit', async () => {
    setApiTokenProvider(() => Promise.resolve(makeJwt({ sub: 'user-only' })))

    await expect(
      createZpeJob({
        queueName: 'default',
        managementNodeId: 'mgmt-node-1',
      }),
    ).rejects.toThrow('Missing required tenant/user claims for runtime submit')
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
