import { FilePanel } from './FilePanel'
import { emptyWorkspaceFile, populatedWorkspaceFile } from './storybook.mocks'

import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'Features/EditorV2/FilePanel',
  component: FilePanel,
  tags: ['autodocs'],
} satisfies Meta<typeof FilePanel>

export default meta

type Story = StoryObj<typeof meta>

export const Default: Story = {
  args: {
    fileId: emptyWorkspaceFile.id,
    data: emptyWorkspaceFile,
    showHeader: true,
  },
  render: (args) => (
    <div className="h-[760px] w-full max-w-6xl">
      <FilePanel {...args} />
    </div>
  ),
}

export const WithData: Story = {
  args: {
    fileId: populatedWorkspaceFile.id,
    data: populatedWorkspaceFile,
    showHeader: true,
  },
  render: (args) => (
    <div className="h-[760px] w-full max-w-6xl">
      <FilePanel {...args} />
    </div>
  ),
}
