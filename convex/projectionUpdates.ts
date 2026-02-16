import { mutation } from './_generated/server'
import { ConvexError, v } from 'convex/values'
import {
  INVALID_PROJECTION_TRANSITION,
  ProjectionTransitionError,
  assertMonotonicProjectionTransition,
  getProjectionStateOrder,
} from './projectionStateMachine'

const projectionStateValidator = v.union(
  v.literal('queued'),
  v.literal('running'),
  v.literal('succeeded'),
  v.literal('failed'),
)

export const applyProjectionUpdate = mutation({
  args: {
    projection_event_id: v.string(),
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
  },
  handler: async (ctx, args) => {
    const target = {
      tenant_id: args.tenant_id,
      workspace_id: args.workspace_id,
      job_id: args.job_id,
      submission_id: args.submission_id,
    }

    const existingReceipt = await (ctx.db
      .query('projection_event_receipts') as any)
      .withIndex('by_projection_event_id', (q: any) =>
        q.eq('projection_event_id', args.projection_event_id),
      )
      .first()

    if (existingReceipt) {
      const isSameTarget =
        existingReceipt.tenant_id === target.tenant_id &&
        existingReceipt.workspace_id === target.workspace_id &&
        existingReceipt.job_id === target.job_id &&
        existingReceipt.submission_id === target.submission_id
      const isSameState =
        existingReceipt.projection_state === args.projection_state

      if (!isSameTarget || !isSameState) {
        throw new ConvexError({
          code: 'projection_event_id_conflict',
          status: 409,
          message: 'projection_event_id was already recorded with different payload',
          existing_target: {
            tenant_id: existingReceipt.tenant_id,
            workspace_id: existingReceipt.workspace_id,
            job_id: existingReceipt.job_id,
            submission_id: existingReceipt.submission_id,
          },
          incoming_target: target,
          existing_projection_state: existingReceipt.projection_state,
          incoming_projection_state: args.projection_state,
        })
      }

      return {
        ok: true,
        applied: false,
        idempotent: true,
        projection_state: existingReceipt.projection_state,
      } as const
    }

    const existingProjection = await (ctx.db.query('projection_states') as any)
      .withIndex('by_target', (q: any) =>
        q
          .eq('tenant_id', args.tenant_id)
          .eq('workspace_id', args.workspace_id)
          .eq('job_id', args.job_id)
          .eq('submission_id', args.submission_id),
      )
      .first()

    const previousState = existingProjection?.projection_state ?? null
    const nowIso = new Date().toISOString()

    if (previousState != null) {
      const previousOrder = getProjectionStateOrder(previousState)
      const nextOrder = getProjectionStateOrder(args.projection_state)

      if (previousState === args.projection_state || nextOrder < previousOrder) {
        await ctx.db.insert('projection_event_receipts', {
          projection_event_id: args.projection_event_id,
          tenant_id: args.tenant_id,
          workspace_id: args.workspace_id,
          job_id: args.job_id,
          submission_id: args.submission_id,
          projection_state: args.projection_state,
          recorded_at: nowIso,
        })
        return {
          ok: true,
          applied: false,
          idempotent: true,
          projection_state: previousState,
        } as const
      }
    }

    try {
      assertMonotonicProjectionTransition(previousState, args.projection_state)
    } catch (error) {
      if (error instanceof ProjectionTransitionError) {
        throw new ConvexError({
          code: INVALID_PROJECTION_TRANSITION,
          status: 409,
          message: error.message,
          previous_projection_state: error.previous_state,
          next_projection_state: error.next_state,
        })
      }
      throw error
    }

    const projectionDocument = {
      tenant_id: args.tenant_id,
      workspace_id: args.workspace_id,
      job_id: args.job_id,
      submission_id: args.submission_id,
      projection_state: args.projection_state,
      source_execution_state: args.source_execution_state,
      occurred_at: args.occurred_at,
      trace_id: args.trace_id,
      display_message: args.display_message,
      result_ref: args.result_ref,
      failure: args.failure,
      last_projection_event_id: args.projection_event_id,
      updated_at: nowIso,
    }

    if (existingProjection) {
      await ctx.db.patch(existingProjection._id, projectionDocument)
    } else {
      await ctx.db.insert('projection_states', projectionDocument)
    }

    await ctx.db.insert('projection_event_receipts', {
      projection_event_id: args.projection_event_id,
      tenant_id: args.tenant_id,
      workspace_id: args.workspace_id,
      job_id: args.job_id,
      submission_id: args.submission_id,
      projection_state: args.projection_state,
      recorded_at: nowIso,
    })

    return {
      ok: true,
      applied: true,
      idempotent: false,
      projection_state: args.projection_state,
    } as const
  },
})
