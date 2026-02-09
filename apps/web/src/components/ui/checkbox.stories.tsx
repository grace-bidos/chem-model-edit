import { Checkbox } from './checkbox'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Checkbox',
  component: Checkbox,
  tags: ['autodocs'],
} satisfies Meta<typeof Checkbox>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
