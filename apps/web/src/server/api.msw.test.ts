import { HttpResponse, http } from 'msw'
import { describe, expect, it, vi } from 'vitest'

import { requestApiInternal } from './api'
import { server } from '@/test/msw/server'

const toBase64Url = (value: string): string => {
  const base64 = btoa(value)
  return base64.replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/g, '')
}

const unsignedJwt = (payload: Record<string, unknown>): string => {
  const header = toBase64Url(JSON.stringify({ alg: 'none', typ: 'JWT' }))
  const body = toBase64Url(JSON.stringify(payload))
  return `${header}.${body}.`
}

describe('requestApiInternal with MSW', () => {
  it('requests JSON payload and forwards auth header', async () => {
    server.use(
      http.get('http://localhost:8000/api/jobs', ({ request }) => {
        return HttpResponse.json({
          auth: request.headers.get('authorization'),
          tenantId: request.headers.get('x-tenant-id'),
          result: 'ok',
        })
      }),
    )

    const response = await requestApiInternal<{
      auth: string | null
      tenantId: string | null
      result: string
    }>({
      path: '/jobs',
      token: 'token-msw',
    })

    expect(response).toEqual({
      auth: 'Bearer token-msw',
      tenantId: null,
      result: 'ok',
    })
  })

  it('forwards explicit tenant header when tenantId is provided', async () => {
    server.use(
      http.get('http://localhost:8000/api/jobs-tenant-explicit', ({ request }) => {
        return HttpResponse.json({
          tenantId: request.headers.get('x-tenant-id'),
        })
      }),
    )

    const response = await requestApiInternal<{
      tenantId: string | null
    }>({
      path: '/jobs-tenant-explicit',
      tenantId: 'tenant-explicit-1',
    })

    expect(response).toEqual({
      tenantId: 'tenant-explicit-1',
    })
  })

  it('derives tenant header from JWT org_id claim when tenantId is omitted', async () => {
    server.use(
      http.get('http://localhost:8000/api/jobs-tenant-derived', ({ request }) => {
        return HttpResponse.json({
          tenantId: request.headers.get('x-tenant-id'),
        })
      }),
    )

    const token = unsignedJwt({ sub: 'user_1', org_id: 'org_abc123' })
    const response = await requestApiInternal<{
      tenantId: string | null
    }>({
      path: '/jobs-tenant-derived',
      token,
    })

    expect(response).toEqual({
      tenantId: 'org_abc123',
    })
  })

  it('supports text response mode', async () => {
    server.use(
      http.get('http://localhost:8000/api/version', () =>
        HttpResponse.text('2026.02'),
      ),
    )

    const response = await requestApiInternal<string>({
      path: '/version',
      responseType: 'text',
    })

    expect(response).toBe('2026.02')
  })

  it('sends JSON body with content-type when body is present', async () => {
    server.use(
      http.post('http://localhost:8000/api/jobs', async ({ request }) => {
        return HttpResponse.json({
          contentType: request.headers.get('content-type'),
          payload: await request.json(),
        })
      }),
    )

    const response = await requestApiInternal<{
      contentType: string | null
      payload: unknown
    }>({
      path: '/jobs',
      method: 'POST',
      body: { mode: 'scf' },
    })

    expect(response).toEqual({
      contentType: 'application/json',
      payload: { mode: 'scf' },
    })
  })

  it('returns undefined for empty successful JSON response', async () => {
    server.use(
      http.get('http://localhost:8000/api/empty', () =>
        new HttpResponse(null, { status: 200 }),
      ),
    )

    const response = await requestApiInternal<undefined>({
      path: '/empty',
    })

    expect(response).toBeUndefined()
  })

  it('uses API envelope message for JSON error responses', async () => {
    server.use(
      http.get('http://localhost:8000/api/error', () =>
        HttpResponse.json(
          {
            error: {
              message: 'upstream failed',
            },
          },
          { status: 500, statusText: 'Internal Server Error' },
        ),
      ),
    )

    await expect(
      requestApiInternal({
        path: '/error',
      }),
    ).rejects.toThrow('upstream failed')
  })

  it('falls back to status text when error body is not JSON', async () => {
    server.use(
      http.get('http://localhost:8000/api/error-text', () =>
        HttpResponse.text('not-json', { status: 502, statusText: 'Bad Gateway' }),
      ),
    )

    await expect(
      requestApiInternal({
        path: '/error-text',
      }),
    ).rejects.toThrow('Bad Gateway')
  })

  it('handles text mode errors with JSON envelope', async () => {
    server.use(
      http.get('http://localhost:8000/api/error-text-mode', () =>
        HttpResponse.json(
          {
            error: {
              message: 'forbidden',
            },
          },
          { status: 403, statusText: 'Forbidden' },
        ),
      ),
    )

    await expect(
      requestApiInternal({
        path: '/error-text-mode',
        responseType: 'text',
      }),
    ).rejects.toThrow('forbidden')
  })

  it('forwards modal proxy auth headers when configured', async () => {
    vi.stubEnv('MODAL_PROXY_KEY', 'proxy-key-test')
    vi.stubEnv('MODAL_PROXY_SECRET', 'proxy-secret-test')

    server.use(
      http.get('http://localhost:8000/api/proxy-auth', ({ request }) => {
        return HttpResponse.json({
          modalKey: request.headers.get('modal-key'),
          modalSecret: request.headers.get('modal-secret'),
        })
      }),
    )

    const response = await requestApiInternal<{
      modalKey: string | null
      modalSecret: string | null
    }>({
      path: '/proxy-auth',
    })

    expect(response).toEqual({
      modalKey: 'proxy-key-test',
      modalSecret: 'proxy-secret-test',
    })
    vi.unstubAllEnvs()
  })
})
