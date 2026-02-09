import { ToolPanel } from './ToolPanel'
import { mockWorkspaceFiles } from './storybook.mocks'

import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'Features/EditorV2/ToolPanel',
  component: ToolPanel,
  tags: ['autodocs'],
} satisfies Meta<typeof ToolPanel>

export default meta

type Story = StoryObj<typeof meta>

export const Default: Story = {
  args: {
    mode: 'supercell',
    structures: [],
    showHeader: true,
    showClose: false,
  },
  render: (args) => (
    <div className="h-[860px] w-full max-w-6xl">
      <ToolPanel {...args} />
    </div>
  ),
}

export const WithData: Story = {
  args: {
    mode: 'supercell',
    structures: mockWorkspaceFiles,
    showHeader: true,
    showClose: false,
  },
  render: (args) => (
    <div className="h-[860px] w-full max-w-6xl">
      <ToolPanel {...args} />
    </div>
  ),
}
