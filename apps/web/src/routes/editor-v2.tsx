import { Suspense, lazy, useEffect, useState } from 'react'
import { createFileRoute } from '@tanstack/react-router'

const EditorV2Page = lazy(() => import('@/features/editor-v2/EditorV2Page'))

export const Route = createFileRoute('/editor-v2')({
  component: EditorV2Route,
})

function EditorV2Route() {
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
      {mounted ? <EditorV2Page /> : fallback}
    </Suspense>
  )
}
