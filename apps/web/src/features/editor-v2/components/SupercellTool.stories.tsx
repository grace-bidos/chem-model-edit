import { SupercellTool } from './SupercellTool'
import { mockWorkspaceFiles } from './storybook.mocks'

import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'Features/EditorV2/SupercellTool',
  component: SupercellTool,
  tags: ['autodocs'],
} satisfies Meta<typeof SupercellTool>

export default meta

type Story = StoryObj<typeof meta>

export const Default: Story = {
  args: {
    structures: [],
  },
}

export const WithData: Story = {
  args: {
    structures: mockWorkspaceFiles,
  },
}
