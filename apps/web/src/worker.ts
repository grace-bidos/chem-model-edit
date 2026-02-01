import { createStartHandler, defaultStreamHandler } from '@tanstack/react-start/server'

const handler = createStartHandler(defaultStreamHandler)

const isAssetRequest = (request: Request): boolean => {
  const { pathname } = new URL(request.url)
  if (pathname.startsWith('/assets/')) {
    return true
  }
  if (
    pathname === '/favicon.ico' ||
    pathname === '/robots.txt' ||
    pathname === '/sitemap.xml'
  ) {
    return true
  }
  return false
}

type AssetsBinding = {
  fetch: (request: Request) => Promise<Response>
}

export default {
  async fetch(request: Request, env: { ASSETS?: AssetsBinding }): Promise<Response> {
    if ((request.method === 'GET' || request.method === 'HEAD') && env.ASSETS) {
      if (isAssetRequest(request)) {
        const assetResponse = await env.ASSETS.fetch(request)
        if (assetResponse.status !== 404) {
          return assetResponse
        }
      }
    }

    return handler(request)
  },
}
