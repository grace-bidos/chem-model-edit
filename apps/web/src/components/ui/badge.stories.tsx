import { Badge } from './badge'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Badge',
  component: Badge,
  tags: ['autodocs'],
} satisfies Meta<typeof Badge>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
