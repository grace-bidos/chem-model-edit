import { createServerFn } from '@tanstack/react-start'

/** API プロキシ呼び出し時に渡すリクエスト情報。 */
export type ApiRequest = {
  path: string
  method?: 'GET' | 'POST' | 'PUT' | 'PATCH' | 'DELETE'
  body?: unknown
  token?: string | null
  tenantId?: string | null
  responseType?: 'json' | 'text'
}

type ApiError = {
  error?: {
    code?: string
    message?: string
    details?: unknown
  }
}

/** API ベースURLを `/api` 終端へ正規化する。 */
const normalizeApiBase = (base: string) => {
  const trimmed = base.endsWith('/') ? base.slice(0, -1) : base
  return trimmed.endsWith('/api') ? trimmed : `${trimmed}/api`
}

const API_BASE = normalizeApiBase(
  process.env.API_BASE_PUBLIC ??
    process.env.API_BASE ??
    import.meta.env.VITE_API_BASE ??
    'http://localhost:8000',
)

const resolveModalProxyHeaders = (): { key: string; secret: string } | null => {
  const key = process.env.MODAL_PROXY_KEY?.trim() ?? ''
  const secret = process.env.MODAL_PROXY_SECRET?.trim() ?? ''
  if (!key || !secret) {
    return null
  }
  return { key, secret }
}

const TENANT_ID_PATTERN = /^[A-Za-z0-9][A-Za-z0-9._:-]{1,127}$/

const normalizeTenantId = (value: unknown): string | null => {
  if (typeof value !== 'string') {
    return null
  }
  const trimmed = value.trim()
  return TENANT_ID_PATTERN.test(trimmed) ? trimmed : null
}

const decodeJwtPayload = (token: string): Record<string, unknown> | null => {
  const parts = token.split('.')
  if (parts.length < 2) {
    return null
  }
  try {
    const payloadPart = parts[1] ?? ''
    const normalized = payloadPart.replace(/-/g, '+').replace(/_/g, '/')
    const padded = normalized.padEnd(
      normalized.length + ((4 - (normalized.length % 4)) % 4),
      '=',
    )
    const payload = JSON.parse(
      atob(padded),
    ) as unknown
    return payload && typeof payload === 'object'
      ? (payload as Record<string, unknown>)
      : null
  } catch (_error) {
    return null
  }
}

const resolveTenantId = (data: ApiRequest): string | null => {
  const explicit = normalizeTenantId(data.tenantId)
  if (explicit) {
    return explicit
  }
  if (!data.token) {
    return null
  }
  const claims = decodeJwtPayload(data.token)
  if (!claims) {
    return null
  }
  const direct =
    normalizeTenantId(claims.tenant_id) ?? normalizeTenantId(claims.tenantId)
  if (direct) {
    return direct
  }
  return normalizeTenantId(claims.org_id) ?? normalizeTenantId(claims.orgId)
}

async function handleResponse<T>(response: Response): Promise<T> {
  if (!response.ok) {
    let message = response.statusText
    try {
      const data = (await response.json()) as ApiError
      if (data.error?.message) {
        message = data.error.message
      }
    } catch (_err) {
      // ignore JSON parsing errors
    }
    throw new Error(message)
  }
  const text = await response.text()
  if (!text) {
    return undefined as T
  }
  return JSON.parse(text) as T
}

async function handleTextResponse(response: Response): Promise<string> {
  if (!response.ok) {
    let message = response.statusText
    try {
      const data = (await response.json()) as ApiError
      if (data.error?.message) {
        message = data.error.message
      }
    } catch (_err) {
      // ignore JSON parsing errors
    }
    throw new Error(message)
  }
  return response.text()
}

type RequestApiFn = ((options: { data: ApiRequest }) => Promise<unknown>) & {
  url: string
}

export async function requestApiInternal<T = unknown>(
  data: ApiRequest,
): Promise<T> {
  const { path, method, body, token } = data
  const responseType = data.responseType ?? 'json'
  const resolvedMethod = method ?? 'GET'

  const headers: Record<string, string> = {}
  if (body !== undefined) {
    headers['Content-Type'] = 'application/json'
  }
  if (token) {
    headers.Authorization = `Bearer ${token}`
  }
  const tenantId = resolveTenantId(data)
  if (tenantId) {
    headers['x-tenant-id'] = tenantId
  }
  const modalProxyHeaders = resolveModalProxyHeaders()
  if (modalProxyHeaders) {
    headers['Modal-Key'] = modalProxyHeaders.key
    headers['Modal-Secret'] = modalProxyHeaders.secret
  }

  const response = await fetch(`${API_BASE}${path}`, {
    method: resolvedMethod,
    headers,
    body: body !== undefined ? JSON.stringify(body) : undefined,
  })

  if (responseType === 'text') {
    return (await handleTextResponse(response)) as T
  }

  return handleResponse<T>(response)
}

/** サーバー側から API へプロキシする共通リクエスト関数。 */
export const requestApi: RequestApiFn = createServerFn({
  method: 'POST',
})
  .inputValidator((data: ApiRequest) => data)
  .handler(async ({ data }) => {
    return requestApiInternal(data)
  })
