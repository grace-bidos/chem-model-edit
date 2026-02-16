import { describe, expect, it } from 'vitest'
import {
  INVALID_PROJECTION_TRANSITION,
  ProjectionTransitionError,
  applyProjectionUpdate,
  createEmptyProjectionAggregate,
} from './projectionStateMachine'

describe('ProjectionUpdate monotonicity and dedup', () => {
  it('applies monotonic transitions in order', () => {
    const aggregate = createEmptyProjectionAggregate()

    const queued = applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-1',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'queued',
    })

    const running = applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-2',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'running',
    })

    const succeeded = applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-3',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'succeeded',
    })

    expect(queued).toEqual({
      applied: true,
      idempotent: false,
      projection_state: 'queued',
    })
    expect(running).toEqual({
      applied: true,
      idempotent: false,
      projection_state: 'running',
    })
    expect(succeeded).toEqual({
      applied: true,
      idempotent: false,
      projection_state: 'succeeded',
    })
  })

  it('treats duplicate projection_event_id as idempotent no-op', () => {
    const aggregate = createEmptyProjectionAggregate()

    const first = applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-replay',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'running',
    })

    const replay = applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-replay',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'failed',
    })

    expect(first).toEqual({
      applied: true,
      idempotent: false,
      projection_state: 'running',
    })
    expect(replay).toEqual({
      applied: false,
      idempotent: true,
      projection_state: 'running',
    })
  })

  it('rejects non-monotonic transitions with stable error semantics', () => {
    const aggregate = createEmptyProjectionAggregate()

    applyProjectionUpdate(aggregate, {
      projection_event_id: 'evt-10',
      tenant_id: 'tenant-1',
      workspace_id: 'ws-1',
      job_id: 'job-1',
      submission_id: 'sub-1',
      projection_state: 'succeeded',
    })

    let thrown
    try {
      applyProjectionUpdate(aggregate, {
        projection_event_id: 'evt-11',
        tenant_id: 'tenant-1',
        workspace_id: 'ws-1',
        job_id: 'job-1',
        submission_id: 'sub-1',
        projection_state: 'running',
      })
    } catch (error) {
      thrown = error
    }

    expect(thrown).toBeInstanceOf(ProjectionTransitionError)
    expect(thrown.code).toBe(INVALID_PROJECTION_TRANSITION)
    expect(thrown.status).toBe(409)
    expect(thrown.previous_state).toBe('succeeded')
    expect(thrown.next_state).toBe('running')
    expect(thrown.message).toBe('invalid projection transition: succeeded -> running')
  })
})
