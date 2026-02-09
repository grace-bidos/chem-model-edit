import { Table } from './table'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'UI/Table',
  component: Table,
  tags: ['autodocs'],
} satisfies Meta<typeof Table>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
