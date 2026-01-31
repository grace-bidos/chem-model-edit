import { Suspense, lazy } from 'react'
import { HeadContent, Outlet, Scripts, createRootRoute } from '@tanstack/react-router'

import Header from '../components/Header'
import '../styles.css'

import type { ReactNode } from 'react'

const runtimeApiBase = import.meta.env.SSR
  ? process.env.API_BASE_PUBLIC ??
    process.env.API_BASE ??
    import.meta.env.VITE_API_BASE ??
    'http://localhost:8000'
  : null

const RouterDevtoolsPanel = import.meta.env.DEV
  ? lazy(() =>
      import('@tanstack/react-router-devtools').then((module) => ({
        default: module.TanStackRouterDevtoolsPanel,
      })),
    )
  : null

const Devtools = import.meta.env.DEV
  ? lazy(() =>
      import('@tanstack/react-devtools').then((module) => ({
        default: module.TanStackDevtools,
      })),
    )
  : null

export const Route = createRootRoute({
  component: RootLayout,
})

function RootLayout() {
  return (
    <RootDocument>
      <Header />
      <Outlet />
      {Devtools && RouterDevtoolsPanel ? (
        <Suspense fallback={null}>
          <Devtools
            config={{
              position: 'bottom-right',
            }}
            plugins={[
              {
                name: 'Tanstack Router',
                render: <RouterDevtoolsPanel />,
              },
            ]}
          />
        </Suspense>
      ) : null}
    </RootDocument>
  )
}

function RootDocument({ children }: { children: ReactNode }) {
  return (
    <html lang="en">
      <head>
        <meta charSet="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        {runtimeApiBase ? (
          <script
            dangerouslySetInnerHTML={{
              __html: `window.__API_BASE__ = ${JSON.stringify(runtimeApiBase)};`,
            }}
          />
        ) : null}
        <HeadContent />
      </head>
      <body>
        {children}
        <Scripts />
      </body>
    </html>
  )
}
