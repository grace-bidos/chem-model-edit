import { Suspense, lazy } from 'react'
import { ClerkProvider } from '@clerk/clerk-react'
import {
  HeadContent,
  Outlet,
  Scripts,
  createRootRoute,
} from '@tanstack/react-router'

import Header from '../components/Header'
import { ClerkTokenBridge } from '../lib/clerk-token-bridge'
import '../styles.css'

import type { ReactNode } from 'react'

const runtimeApiBase = import.meta.env.SSR
  ? (process.env.API_BASE_PUBLIC ??
    process.env.API_BASE ??
    import.meta.env.VITE_API_BASE ??
    'http://localhost:8000')
  : null

const clerkPublishableKey = import.meta.env.VITE_CLERK_PUBLISHABLE_KEY
const clerkKey = clerkPublishableKey ?? 'pk_test_Y2xlcmsuZXhhbXBsZS5jb20k'

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

/** ルートレイアウトのエントリーポイント。 */
export const Route = createRootRoute({
  component: RootLayout,
})

function RootLayout() {
  const appShell = (
    <>
      <ClerkTokenBridge />
      {!clerkPublishableKey ? (
        <div className="border-b border-amber-300 bg-amber-50 px-3 py-2 text-xs text-amber-800">
          VITE_CLERK_PUBLISHABLE_KEY is not set. Clerk auth will not be usable.
        </div>
      ) : null}
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
    </>
  )

  return (
    <RootDocument>
      <ClerkProvider publishableKey={clerkKey}>{appShell}</ClerkProvider>
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
