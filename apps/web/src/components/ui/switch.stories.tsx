import { Switch } from './switch'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Switch',
  component: Switch,
  tags: ['autodocs'],
} satisfies Meta<typeof Switch>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
