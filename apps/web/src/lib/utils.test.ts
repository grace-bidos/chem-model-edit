import { describe, expect, it } from 'vitest'

import { cn } from './utils'

describe('cn', () => {
  it('merges tailwind utility conflicts', () => {
    expect(cn('p-2', 'p-4', 'text-slate-500', 'text-slate-700')).toBe(
      'p-4 text-slate-700',
    )
  })
})
