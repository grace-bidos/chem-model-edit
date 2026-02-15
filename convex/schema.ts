import { defineSchema, defineTable } from 'convex/server'
import { v } from 'convex/values'

const projectionStateValidator = v.union(
  v.literal('queued'),
  v.literal('running'),
  v.literal('succeeded'),
  v.literal('failed'),
)

export default defineSchema({
  projection_states: defineTable({
    tenant_id: v.string(),
    workspace_id: v.string(),
    job_id: v.string(),
    submission_id: v.string(),
    projection_state: projectionStateValidator,
    source_execution_state: v.union(
      v.literal('accepted'),
      v.literal('running'),
      v.literal('completed'),
      v.literal('failed'),
    ),
    occurred_at: v.string(),
    trace_id: v.string(),
    display_message: v.optional(v.string()),
    result_ref: v.optional(
      v.object({
        output_uri: v.string(),
        summary: v.optional(v.string()),
      }),
    ),
    failure: v.optional(
      v.object({
        code: v.string(),
        message: v.string(),
      }),
    ),
    last_projection_event_id: v.string(),
    updated_at: v.string(),
  }).index('by_target', [
    'tenant_id',
    'workspace_id',
    'job_id',
    'submission_id',
  ]),

  projection_event_receipts: defineTable({
    projection_event_id: v.string(),
    tenant_id: v.string(),
    workspace_id: v.string(),
    job_id: v.string(),
    submission_id: v.string(),
    projection_state: projectionStateValidator,
    recorded_at: v.string(),
  }).index('by_projection_event_id', ['projection_event_id']),
})
