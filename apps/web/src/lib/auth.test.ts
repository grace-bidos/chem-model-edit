/* @vitest-environment jsdom */
import { describe, expect, it, vi } from 'vitest'

import {
  clearSession,
  getAuthToken,
  getStoredSession,
  storeSession,
} from './auth'

describe('auth session storage', () => {
  it('stores and reads session from localStorage', () => {
    clearSession()
    storeSession({
      token: 'token-1',
      expires_at: '2099-01-01T00:00:00Z',
      user: {
        id: 'user-1',
        email: 'user-1@example.com',
        created_at: '2026-01-01T00:00:00Z',
      },
    })

    expect(getStoredSession()).toEqual({
      token: 'token-1',
      expires_at: '2099-01-01T00:00:00Z',
      user: {
        id: 'user-1',
        email: 'user-1@example.com',
        created_at: '2026-01-01T00:00:00Z',
      },
    })
    expect(getAuthToken()).toBe('token-1')
  })

  it('returns null for invalid JSON payload', () => {
    window.localStorage.setItem('chem-model-auth-session', '{')
    expect(getStoredSession()).toBeNull()
  })

  it('returns null/no-op when window is unavailable', () => {
    vi.stubGlobal('window', undefined)

    try {
      expect(getStoredSession()).toBeNull()
      expect(getAuthToken()).toBeNull()
      expect(() =>
        storeSession({
          token: 'token-2',
          expires_at: '2099-01-01T00:00:00Z',
          user: {
            id: 'user-2',
            email: 'user-2@example.com',
            created_at: '2026-01-01T00:00:00Z',
          },
        }),
      ).not.toThrow()
      expect(() => clearSession()).not.toThrow()
    } finally {
      vi.unstubAllGlobals()
    }
  })
})
