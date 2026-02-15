export const JOB_STATES = ['queued', 'started', 'finished', 'failed'] as const

export type JobState = (typeof JOB_STATES)[number]

export const TERMINAL_JOB_STATES = ['finished', 'failed'] as const

export const canTransitionJobState = (
  previous: JobState | null | undefined,
  next: JobState,
): boolean => {
  if (previous == null) return true
  if (previous === next) return true

  switch (previous) {
    case 'finished':
      return false
    case 'failed':
      return next === 'queued'
    case 'queued':
      return next === 'started' || next === 'failed' || next === 'finished'
    case 'started':
      return next === 'queued' || next === 'finished' || next === 'failed'
    default:
      return false
  }
}
