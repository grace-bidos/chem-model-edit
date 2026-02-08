import { createStart } from '@tanstack/react-start'

/** SPA モードで起動する TanStack Start インスタンス。 */
export const startInstance = createStart(() => ({
  defaultSsr: false,
}))
