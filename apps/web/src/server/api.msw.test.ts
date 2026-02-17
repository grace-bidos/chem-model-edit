import { HttpResponse, http } from 'msw'
import { describe, expect, it } from 'vitest'

import { requestApiInternal } from './api'
import { server } from '@/test/msw/server'

describe('requestApiInternal with MSW', () => {
  it('requests JSON payload and forwards auth header', async () => {
    server.use(
      http.get('http://localhost:8000/api/jobs', ({ request }) => {
        return HttpResponse.json({
          auth: request.headers.get('authorization'),
          result: 'ok',
        })
      }),
    )

    const response = await requestApiInternal<{
      auth: string | null
      result: string
    }>({
      path: '/jobs',
      token: 'token-msw',
    })

    expect(response).toEqual({
      auth: 'Bearer token-msw',
      result: 'ok',
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
})
