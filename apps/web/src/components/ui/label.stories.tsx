import { Label } from './label'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Label',
  component: Label,
  tags: ['autodocs'],
} satisfies Meta<typeof Label>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
