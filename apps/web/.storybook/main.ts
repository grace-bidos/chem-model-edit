import { mergeConfig } from 'vite'
import tailwindcss from '@tailwindcss/vite'
import viteTsConfigPaths from 'vite-tsconfig-paths'
import type { StorybookConfig } from '@storybook/react-vite'

const config: StorybookConfig = {
  stories: ['../src/**/*.stories.@(js|jsx|ts|tsx|mdx)'],
  addons: [
    '@storybook/addon-links',
    '@storybook/addon-essentials',
    '@storybook/addon-interactions',
  ],
  framework: {
    name: '@storybook/react-vite',
    options: {},
  },
  docs: {
    autodocs: 'tag',
  },
  viteFinal: (viteConfig) =>
    mergeConfig(viteConfig, {
      plugins: [
        viteTsConfigPaths({
          projects: ['../tsconfig.json'],
        }),
        tailwindcss(),
      ],
    }),
}

export default config
