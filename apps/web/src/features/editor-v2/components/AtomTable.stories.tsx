import { AtomTable } from './AtomTable'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'Features/EditorV2/AtomTable',
  component: AtomTable,
  tags: ['autodocs'],
} satisfies Meta<typeof AtomTable>

export default meta

type Story = StoryObj<typeof meta>

const rows = [
  { index: 0, symbol: 'Si', x: 0, y: 0, z: 0 },
  { index: 1, symbol: 'O', x: 1.6, y: 0, z: 0 },
  { index: 2, symbol: 'O', x: -1.6, y: 0, z: 0 },
]

export const Default: Story = {
  args: {
    rows,
    selectionEnabled: false,
  },
}

export const WithData: Story = {
  args: {
    rows,
    selectedIndices: new Set([1]),
    fixedIndices: new Set([2]),
    selectionEnabled: true,
    onRowClick: () => {},
  },
}
