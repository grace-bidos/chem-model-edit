import { CollapsibleSection } from './CollapsibleSection'
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: 'Features/EditorV2/CollapsibleSection',
  component: CollapsibleSection,
  tags: ['autodocs'],
} satisfies Meta<typeof CollapsibleSection>

export default meta

type Story = StoryObj<typeof meta>

export const Default: Story = {
  args: {
    title: 'Supercell Options',
    defaultOpen: false,
    children: 'Section content goes here.',
  },
  render: (args) => (
    <CollapsibleSection {...args}>
      <p className="text-sm text-slate-600">Section content goes here.</p>
    </CollapsibleSection>
  ),
}

export const WithData: Story = {
  args: {
    title: 'Selected Structures',
    defaultOpen: true,
    children: 'Selected structures list',
  },
  render: (args) => (
    <CollapsibleSection {...args}>
      <ul className="space-y-1 text-sm text-slate-700">
        <li>SiO slab (3 atoms)</li>
        <li>CH4 molecule (5 atoms)</li>
      </ul>
    </CollapsibleSection>
  ),
}
