import { Navigate, createFileRoute } from '@tanstack/react-router'

/** ルート `/` をエディタ画面へリダイレクトする。 */
export const Route = createFileRoute('/')({ component: Index })

function Index() {
  return <Navigate to="/editor" />
}
