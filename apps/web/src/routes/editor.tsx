import { Suspense, lazy, useEffect, useState } from 'react'
import { createFileRoute } from '@tanstack/react-router'

const EditorPage = lazy(() => import('@/features/editor-v2/EditorV2Page'))

/** エディタ画面ルート。 */
export const Route = createFileRoute('/editor')({
  component: EditorRoute,
})

function EditorRoute() {
  const [mounted, setMounted] = useState(false)

  useEffect(() => {
    setMounted(true)
  }, [])

  const fallback = (
    <div className="flex min-h-[40vh] items-center justify-center text-sm text-slate-500">
      Loading editor...
    </div>
  )

  return (
    <Suspense fallback={fallback}>
      {mounted ? <EditorPage /> : fallback}
    </Suspense>
  )
}
