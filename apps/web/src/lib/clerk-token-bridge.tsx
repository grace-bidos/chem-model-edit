import { useAuth } from '@clerk/clerk-react'
import { useEffect } from 'react'

import { setApiTokenProvider } from './api'

/** Clerk のセッショントークンを API クライアントへ橋渡しする。 */
export function ClerkTokenBridge() {
  const { isSignedIn, getToken } = useAuth()

  useEffect(() => {
    setApiTokenProvider(async () => {
      if (!isSignedIn) {
        return null
      }
      return getToken()
    })
    return () => {
      setApiTokenProvider(null)
    }
  }, [getToken, isSignedIn])

  return null
}
