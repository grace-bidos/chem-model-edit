import { Slider } from './slider'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Slider',
  component: Slider,
  tags: ['autodocs'],
} satisfies Meta<typeof Slider>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
