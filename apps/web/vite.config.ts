import { defineConfig, type PluginOption } from 'vite'
import { devtools } from '@tanstack/devtools-vite'
import { tanstackStart } from '@tanstack/react-start/plugin/vite'
import viteReact from '@vitejs/plugin-react'
import viteTsConfigPaths from 'vite-tsconfig-paths'
import tailwindcss from '@tailwindcss/vite'
import { nitro } from 'nitro/vite'

const config = defineConfig(({ mode }) => {
  const isTest = mode === 'test' || process.env.VITEST === 'true'
  const plugins = [
    !isTest &&
      devtools({
        eventBusConfig: {
          enabled: process.env.TANSTACK_DEVTOOLS_EVENTBUS === 'true',
          port: Number(process.env.TANSTACK_DEVTOOLS_PORT ?? 42069),
        },
      }),
    !isTest && nitro(),
    // this is the plugin that enables path aliases
    viteTsConfigPaths({
      projects: ['./tsconfig.json'],
    }),
    tailwindcss(),
    !isTest && tanstackStart(),
    viteReact(),
  ].filter(Boolean) as PluginOption[]

  return { plugins }
})

export default config
