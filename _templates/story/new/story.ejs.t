---
to: <%= dir %>/<%= fileBase %>.stories.tsx
---
<% if (exportStyle === 'default') { -%>
import <%= componentName %> from '<%= importPath %>'
<% } else { -%>
import { <%= exportName %> } from '<%= importPath %>'
<% } -%>
import type { Meta, StoryObj } from '@storybook/react'

const meta = {
  title: '<%= storyTitle %>',
  component: <%= exportStyle === 'default' ? componentName : exportName %>,
  tags: ['autodocs'],
} satisfies Meta<typeof <%= exportStyle === 'default' ? componentName : exportName %>>

export default meta

type Story = StoryObj<typeof meta>

export const Primary: Story = {}
