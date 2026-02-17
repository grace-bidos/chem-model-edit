import { URL, fileURLToPath } from 'node:url'
import { cloudflare } from '@cloudflare/vite-plugin'
import { devtools } from '@tanstack/devtools-vite'
import { tanstackStart } from '@tanstack/react-start/plugin/vite'
import viteReact from '@vitejs/plugin-react'
import viteTsConfigPaths from 'vite-tsconfig-paths'
import tailwindcss from '@tailwindcss/vite'
import { defineConfig } from 'vite'
import { configDefaults } from 'vitest/config'
import type { PluginOption } from 'vite'

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
    // this is the plugin that enables path aliases
    viteTsConfigPaths({
      projects: ['./tsconfig.json'],
    }),
    tailwindcss(),
    !isTest && cloudflare({ viteEnvironment: { name: 'ssr' } }),
    !isTest && tanstackStart({ spa: { enabled: true } }),
    viteReact(),
  ].filter(Boolean) as Array<PluginOption>

  return {
    plugins,
    resolve: {
      alias: {
        '@': fileURLToPath(new URL('./src', import.meta.url)),
      },
    },
    optimizeDeps: {
      include: ['dockview', 'dockview-react', 'dockview-core'],
    },
    test: {
      setupFiles: ['./src/test/setup.ts'],
      exclude: [...configDefaults.exclude, 'e2e/**', '.stryker-tmp/**'],
      coverage: {
        provider: 'v8',
        reporter: ['text', 'html', 'json-summary'],
        reportsDirectory: './coverage',
        include: [
          'src/lib/auth.ts',
          'src/lib/utils.ts',
          'src/lib/xyz.ts',
          'src/server/api.ts',
          'src/server/zpe-aggregation.ts',
          'src/components/compare/align.ts',
          'src/features/editor-v2/components/AtomTable.tsx',
          'src/features/editor-v2/components/FilePanel.tsx',
        ],
        exclude: [
          'src/**/*.stories.*',
          'src/**/*.test.*',
          'src/**/*.a11y.test.*',
          'src/**/*.fastcheck.test.*',
          'src/routeTree.gen.ts',
        ],
        thresholds: {
          lines: 80,
          functions: 80,
          branches: 70,
          statements: 80,
        },
      },
    },
  }
})

export default config
