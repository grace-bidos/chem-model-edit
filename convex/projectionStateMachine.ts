export const PROJECTION_STATES = ['queued', 'running', 'succeeded', 'failed'] as const

export type ProjectionState = (typeof PROJECTION_STATES)[number]

export interface ProjectionUpdateInput {
  projection_event_id: string
  tenant_id: string
  workspace_id: string
  job_id: string
  submission_id: string
  projection_state: ProjectionState
}

export interface ProjectionRecord {
  target_key: string
  projection_state: ProjectionState
}

export interface ProjectionAggregate {
  projections_by_target: Map<string, ProjectionRecord>
  seen_projection_event_ids: Set<string>
}

export interface ApplyProjectionUpdateResult {
  applied: boolean
  idempotent: boolean
  projection_state: ProjectionState
}

export const INVALID_PROJECTION_TRANSITION =
  'invalid_projection_transition' as const

const PROJECTION_STATE_ORDER: Record<ProjectionState, number> = {
  queued: 0,
  running: 1,
  succeeded: 2,
  failed: 2,
}

const TERMINAL_PROJECTION_STATES = new Set<ProjectionState>([
  'succeeded',
  'failed',
])

export class ProjectionTransitionError extends Error {
  readonly code = INVALID_PROJECTION_TRANSITION
  readonly status = 409

  constructor(
    readonly previous_state: ProjectionState,
    readonly next_state: ProjectionState,
  ) {
    super(
      `invalid projection transition: ${previous_state} -> ${next_state}`,
    )
    this.name = 'ProjectionTransitionError'
  }
}

export const buildProjectionTargetKey = (
  input: Pick<
    ProjectionUpdateInput,
    'tenant_id' | 'workspace_id' | 'job_id' | 'submission_id'
  >,
): string => {
  return [
    input.tenant_id,
    input.workspace_id,
    input.job_id,
    input.submission_id,
  ].join('::')
}

export const isMonotonicProjectionTransition = (
  previous: ProjectionState | null | undefined,
  next: ProjectionState,
): boolean => {
  if (previous == null) return true
  if (previous === next) return true
  if (TERMINAL_PROJECTION_STATES.has(previous)) return false

  return PROJECTION_STATE_ORDER[next] >= PROJECTION_STATE_ORDER[previous]
}

export const assertMonotonicProjectionTransition = (
  previous: ProjectionState | null | undefined,
  next: ProjectionState,
): void => {
  if (previous == null) return
  if (isMonotonicProjectionTransition(previous, next)) return

  throw new ProjectionTransitionError(previous, next)
}

export const createEmptyProjectionAggregate = (): ProjectionAggregate => {
  return {
    projections_by_target: new Map<string, ProjectionRecord>(),
    seen_projection_event_ids: new Set<string>(),
  }
}

export const applyProjectionUpdate = (
  aggregate: ProjectionAggregate,
  update: ProjectionUpdateInput,
): ApplyProjectionUpdateResult => {
  if (aggregate.seen_projection_event_ids.has(update.projection_event_id)) {
    const targetKey = buildProjectionTargetKey(update)
    const current = aggregate.projections_by_target.get(targetKey)
    return {
      applied: false,
      idempotent: true,
      projection_state: current?.projection_state ?? update.projection_state,
    }
  }

  const targetKey = buildProjectionTargetKey(update)
  const current = aggregate.projections_by_target.get(targetKey)
  const previousState = current?.projection_state ?? null

  assertMonotonicProjectionTransition(previousState, update.projection_state)

  aggregate.projections_by_target.set(targetKey, {
    target_key: targetKey,
    projection_state: update.projection_state,
  })
  aggregate.seen_projection_event_ids.add(update.projection_event_id)

  return {
    applied: true,
    idempotent: false,
    projection_state: update.projection_state,
  }
}
