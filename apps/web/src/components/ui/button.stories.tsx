import { Button } from './button'

import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Button',
  component: Button,
  args: {
    children: 'Button',
  },
} satisfies Meta<typeof Button>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}

export const Destructive: Story = {
  args: {
    variant: 'destructive',
  },
}

export const Outline: Story = {
  args: {
    variant: 'outline',
  },
}

export const Small: Story = {
  args: {
    size: 'sm',
  },
}
