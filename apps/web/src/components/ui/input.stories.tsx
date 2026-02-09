import { Input } from './input'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Input',
  component: Input,
  tags: ['autodocs'],
} satisfies Meta<typeof Input>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
