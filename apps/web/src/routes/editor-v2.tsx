import { Suspense, lazy } from 'react'
import { createFileRoute } from '@tanstack/react-router'

const EditorV2Page = lazy(() => import('@/features/editor-v2/EditorV2Page'))

export const Route = createFileRoute('/editor-v2')({
  component: EditorV2Route,
})

function EditorV2Route() {
  if (typeof document === 'undefined') {
    return null
  }
  return (
    <Suspense
      fallback={
        <div className="flex min-h-[40vh] items-center justify-center text-sm text-slate-500">
          Loading editor...
        </div>
      }
    >
      <EditorV2Page />
    </Suspense>
  )
}
