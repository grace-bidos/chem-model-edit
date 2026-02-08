import { createRouter } from '@tanstack/react-router'

import { routeTree } from './routeTree.gen'

/** アプリケーション全体で共有する TanStack Router インスタンスを生成する。 */
export const getRouter = () => {
  const router = createRouter({
    routeTree,
    context: {},

    scrollRestoration: true,
    defaultPreloadStaleTime: 0,
  })

  return router
}

declare module '@tanstack/react-router' {
  interface Register {
    router: ReturnType<typeof getRouter>
  }
}
