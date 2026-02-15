import { createServerFn } from '@tanstack/react-start'

/** API プロキシ呼び出し時に渡すリクエスト情報。 */
export type ApiRequest = {
  path: string
  method?: 'GET' | 'POST' | 'PUT' | 'PATCH' | 'DELETE'
  body?: unknown
  token?: string | null
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
