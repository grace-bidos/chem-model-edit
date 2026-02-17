//  @ts-check

import { tanstackConfig } from '@tanstack/eslint-config'
import tsdoc from 'eslint-plugin-tsdoc'

export default [
  {
    ignores: [
      '.output/**',
      'dist/**',
      '.tanstack/**',
      '.stryker-tmp/**',
      'storybook-static/**',
    ],
  },
  ...tanstackConfig,
  {
    files: ['**/*.{ts,tsx,js,jsx}'],
    languageOptions: {
      parserOptions: {
        project: ['./tsconfig.eslint.json'],
        tsconfigRootDir: import.meta.dirname,
      },
    },
  },
  {
    files: ['src/**/*.{ts,tsx}'],
    ignores: [
      'src/routeTree.gen.ts',
      'src/**/*.test.ts',
      'src/**/*.test.tsx',
      'src/**/*.stories.ts',
      'src/**/*.stories.tsx',
    ],
    plugins: {
      tsdoc,
    },
    rules: {
      'tsdoc/syntax': 'error',
    },
  },
]
