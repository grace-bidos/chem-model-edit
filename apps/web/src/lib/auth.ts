export type LegacyAuthUser = {
  id: string
  email: string
  created_at: string
}

export type LegacyAuthSession = {
  token: string
  expires_at: string
  user: LegacyAuthUser
}

const SESSION_KEY = 'chem-model-auth-session'

/** ローカルストレージから認証セッションを取得する。 */
export function getStoredSession(): LegacyAuthSession | null {
  if (typeof window === 'undefined') {
    return null
  }
  const raw = window.localStorage.getItem(SESSION_KEY)
  if (!raw) {
    return null
  }
  try {
    return JSON.parse(raw) as LegacyAuthSession
  } catch (_err) {
    return null
  }
}

/** 認証セッションをローカルストレージへ保存する。 */
export function storeSession(session: LegacyAuthSession) {
  if (typeof window === 'undefined') {
    return
  }
  window.localStorage.setItem(SESSION_KEY, JSON.stringify(session))
}

/** ローカルストレージ上の認証セッションを削除する。 */
export function clearSession() {
  if (typeof window === 'undefined') {
    return
  }
  window.localStorage.removeItem(SESSION_KEY)
}

/** 保存済みセッションから Bearer トークンを取得する。 */
export function getAuthToken(): string | null {
  const session = getStoredSession()
  return session?.token ?? null
}
