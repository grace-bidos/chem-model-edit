import { defineSchema, defineTable } from 'convex/server'
import { v } from 'convex/values'

const projectionStateValidator = v.union(
  v.literal('queued'),
  v.literal('running'),
  v.literal('succeeded'),
  v.literal('failed'),
)

export default defineSchema({
  structures: defineTable({
    tenant_id: v.string(),
    workspace_id: v.string(),
    structure_id: v.string(),
    source: v.string(),
    cif: v.string(),
    structure: v.any(),
    params: v.optional(v.any()),
    created_at: v.string(),
    updated_at: v.string(),
    raw_input_chunk_count: v.number(),
  })
    .index("by_key", ["tenant_id", "workspace_id", "structure_id"])
    .index("by_structure_id", ["structure_id"]),

  structure_raw_chunks: defineTable({
    tenant_id: v.string(),
    workspace_id: v.string(),
    structure_id: v.string(),
    chunk_index: v.number(),
    chunk_text: v.string(),
  }).index("by_structure", [
    "tenant_id",
    "workspace_id",
    "structure_id",
    "chunk_index",
  ]),

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
    incoming_projection_state: v.optional(projectionStateValidator),
    recorded_at: v.string(),
  })
    .index('by_projection_event_id', ['projection_event_id'])
    .index('by_target_projection_event_id', [
      'tenant_id',
      'workspace_id',
      'job_id',
      'submission_id',
      'projection_event_id',
    ]),

  runtime_targets: defineTable({
    tenant_id: v.string(),
    user_id: v.string(),
    target_id: v.string(),
    queue_name: v.string(),
    server_id: v.string(),
    registered_at: v.string(),
    name: v.optional(v.string()),
    metadata: v.optional(v.any()),
  })
    .index('by_target_id', ['target_id'])
    .index('by_user', ['tenant_id', 'user_id']),

  runtime_active_targets: defineTable({
    tenant_id: v.string(),
    user_id: v.string(),
    target_id: v.string(),
    updated_at: v.string(),
  }).index('by_user', ['tenant_id', 'user_id']),

  runtime_join_tokens: defineTable({
    token: v.string(),
    tenant_id: v.string(),
    owner_user_id: v.string(),
    queue_name: v.string(),
    created_at: v.string(),
    expires_at: v.string(),
    used_at: v.optional(v.string()),
    node_name_hint: v.optional(v.string()),
    label: v.optional(v.string()),
  }).index('by_token', ['token']),
})
